#!/bin/sh

######
# Description: PIDFE is used to detect P-element insertions and estimate their population frequencies from paired-end reads
# 11/25/2019
# Contact: Shuo Zhang (shuozhang23@gmail.com)
# Copyright (c) 2019 Kelleher Lab at the University of Houston
######


## default options, which are able to be overidden by command line arguments

# the directory that contains and only contains paired-end reads, i.e., *_R1.fastq and *_R2.fastq. Asterisk(*) can represent zero or any number of characters.
inDir=$PWD

# output file, default: p_insertions.txt
output="p_insertions.txt"

usage() {
    echo "Usage: $0 [-i inDir] [-o output] <-g refGenome> <-p refP>"
    echo "Required arguments:"
    echo "  <-g refGenome> defines the reference genome where P-elements reside."
    echo "  <-p refP> defines the consensus P-element sequence."
    echo ""
    echo "Optional arguments:"
    echo "  [-i inDir] is the directory containing paired-end read, i.e., *_R1.fastq and *_R2.fastq. Asterisk(*) can represent zero or any number of characters. Default: current working directoary"
    echo "  [-o output] is the output file name. Default: p_insertions.txt"
    exit $1
}



while getopts "i:o:g:p:h" opt
do
	case $opt in
		i) inDir=$OPTARG ;;
		o) output=$OPTARG ;;
		g) refGenome=$OPTARG ;;
		p) refP=$OPTARG ;;
        h) usage 0;;
        [?]) usage 1;;
	esac
done


### check if required variables are provided
# check if reference genome is provided
if [ ! -e $refGenome ]
then
    echo "*** error: reference genome must be provide or provided reference genome doesn't exist ***"
    usage 1
fi

# check if P-element reference is provided
if [ ! -e $refP ]
then
    echo "*** error: P-element reference genome must be provide or provided P-element reference genome doesn't exist ***"
    usage 1
fi

cd $inDir
read1=*_R1.fastq
read2=*_R2.fastq
# check if read1 and read2 exist
if [ ! -e $read1 ]
then
    echo "*** error: *_R1.fastq doesn't exist in ${inDir} ***"
    usage 1
fi
if [ ! -e $read2 ]
then
    echo "*** error: *_R2.fastq doesn't exist in ${inDir} ***"
    usage 1
fi

echo $PWD
echo `date`


# align paired-end reads to P-element, respectively.
[ ! -d p_alignment ] && mkdir p_alignment
echo "align "$read1" and $read2 to $refP"
echo "cmd: bowtie2 --local -q --no-unal -x $refP -U ${read1} -S p_alignment/read1_p.sam"
echo "cmd: bowtie2 --local -q --no-unal -x $refP -U ${read2} -S p_alignment/read2_p.sam"
bowtie2 --local -q --no-unal -x $refP -U $read1 -S p_alignment/read1_p.sam
bowtie2 --local -q --no-unal -x $refP -U $read2 -S p_alignment/read2_p.sam
echo "Aligning to reference genome is done\n"

# convert .sam to .bed
echo "convert .sam to .bed"
cd p_alignment

###added by erin to create files where reverse complemented alignments are flagged s 'r' I think the original code relied on a function of samtools that no longer exists.
samtools view -h read1_p.sam | perl -pe 's/\t0\tP/\t\tP/g;' |  perl -pe 's/\t0\tP/\tr\tP/g;' >formatted_read1_p.sam
samtools view -h read2_p.sam | perl -pe 's/\t0\tP/\t\tP/g;' |  perl -pe 's/\t0\tP/\tr\tP/g;' >formatted_read2_p.sam

samtools view -bS read1_p.sam > read1_p.bam
bedtools bamtobed -i read1_p.bam -tag NM > read1_p.bed
samtools view -bS read2_p.sam > read2_p.bam
bedtools bamtobed -i read2_p.bam -tag NM > read2_p.bed
rm read1_p.bam read2_p.bam
cd ..
echo "conversion is done"


# pick reads mapped to P-elements and combine read1 and read2 (aligned to P-element) to a table
pick_p_mapped_reads.pl p_alignment/read1_p.bed p_alignment/read2_p.bed p_mapped_read.list
combine_paired_bed.pl p_alignment/read1_p.bed p_alignment/read2_p.bed p_mapped_read.list combined_paired_p.table
rm p_alignment/read1_p.bed p_alignment/read2_p.bed

# pick potential reads
pick_given_reads.pl -list p_mapped_read.list -i $read1 -out potential_read1.fq
pick_given_reads.pl -list p_mapped_read.list -i $read2 -out potential_read2.fq

# pick P-element-trimmed reads
pick_p_trimmed_seq.pl -mapped_list combined_paired_p.table -potential_read1 potential_read1.fq -potential_read2 potential_read2.fq -read1_sam p_alignment/read1_p.sam -read2_sam p_alignment/read2_p.sam -single_read1 single_read1.fq -single_read2 single_read2.fq -paired_read1 paired_read1.fq -paired_read2 paired_read2.fq

# mapping P-derived reads that P-element sequences are trmmmed to reference genome
bowtie2 -q -x $refGenome -1 paired_read1.fq -2 paired_read2.fq -S paired_dm6.sam
bowtie2 -q -x $refGenome -U single_read1.fq,single_read2.fq -S single_dm6.sam

# pick alignments with mapping quality score >= 20 and convert .sam to .bed
samtools view -bSf 0x2 -q 20 paired_dm6.sam > paired_dm6.bam
bedtools bamtobed -ed -i paired_dm6.bam > paired_dm6.bed
samtools view -bS -q 20 single_dm6.sam > single_dm6.bam
bedtools bamtobed -ed -i single_dm6.bam > single_dm6.bed
rm paired_dm6.bam single_dm6.bam


# determine P-element insertion breakpoints
echo "...determining P-element insertion breakpoints starts..."
paired_breakpoint.pl combined_paired_p.table paired_dm6.bed counted_reads.txt paired_insertions.txt
single_breakpoint.pl combined_paired_p.table paired_insertions.txt single_dm6.bed counted_reads.txt paired_single_insertions.txt
combine_insertion_within_line.pl paired_single_insertions.txt combined_paired_single_insertions.txt
recruit_more_paired_reads.pl counted_reads.txt paired_dm6.bed combined_paired_single_insertions.txt combined_paired_p.table recruited_paired_insertions.txt
recruit_more_single_reads.pl counted_reads.txt recruited_paired_insertions.txt combined_paired_p.table single_dm6.bed $output

rm paired_dm6.sam single_dm6.sam paired_dm6.bed single_dm6.bed
rm paired_insertions.txt paired_single_insertions.txt combined_paired_single_insertions.txt recruited_paired_insertions.txt
rm counted_reads.txt combined_paired_p.table p_mapped_read.list
echo "...P-element insertion breakpoints are determined"


# estimate P-element insertion frequency
echo "...estimating insertion frequency begins..."

# build a reference containing constructed pseudo-genomes
echo "...build pseudogenomes..."
mkdir p_pseudo
make_P_pseudo_only.pl -f 500 -g $refGenome -p $refP -i $output -o p_pseudo/p_pseudo.fa
cd p_pseudo
cat p_pseudo.fa $refGenome > combined.fa
rm p_pseudo.fa
mv combined.fa p_pseudo.fa
bowtie2-build -f p_pseudo.fa p_pseudo.fa
cd ..
echo "...pseudo-genomes are built..."

# align paired-end reads to pseudogenome reference, index and sort alignments
echo "...align to pseudo-genomes..."
bowtie2 -p 16 --reorder -q --no-mixed --no-discordant --no-unal -x p_pseudo/p_pseudo.fa -1 $read1 -2 $read2 -S pseudo.sam
echo "...aligning to pseudo-genomes is done..."

###get the header from the samfile
samtools view -h pseudo.sam | grep -P  "^@" >new.pseudo.sam
### get the reads with edit distance <4
samtools view -f 0X2 pseudo.sam | grep -P  "NM:i:[0-3]" >>new.pseudo.sam
####convert to bam
samtools view -b new.pseudo.sam >new.pseudo.bam
####sort reads
samtools sort new.pseudo.bam -o sorted.new.pseudo.bam
#### index reads
samtools index sorted.new.pseudo.bam


### isolate read pairs that do not span breakpoints, convert to fastq
awk '{print $1.":inserted\t1007\t3915 "}' $output >insertions.bed

echo "...collecting reference reads..."
samtools view -@ 8 -U reference.reads.bam -L insertions.bed sorted.new.pseudo.bam >insertion.reads.bam
### sort by read names
samtools sort  -@ 8 -n reference.reads.bam -o sorted.reference.reads.bam
### convert to fastq
samtools bam2fq sorted.reference.reads.bam >reference_reads.fq
paired.pl -f reference_reads.fq


### align read pairs that do not support the insertion to the reference genome
echo "...aligning to reference..."
bowtie2 -p 16 --reorder -q --no-mixed --no-discordant --no-unal -x $refGenome -1 reference_reads_R1.fq -2 reference_reads_R2.fq -S reference.sam
echo "...alignment to reference complete..."
echo "...sorting reference reads..."
samtools view  -@ 8 -h reference.sam | grep -P  "^@" >new.reference.sam
samtools view -@ 8 -f 0X2 reference.sam | grep -P  "NM:i:[0-3]" >>new.reference.sam
samtools view  -@ 8 -b new.reference.sam >new.reference.bam
samtools sort -@ 8 new.reference.bam -o sorted.new.reference.bam
samtools index sorted.new.reference.bam

echo "...calculating frequencies.."
newfreq.pl -f 1000 -i $output

echo "...aligning to pseudo-genomes is done..."




echo "...PIDFE is done..."
echo `date`



