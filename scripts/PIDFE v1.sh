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

<<COMMENT
# align paired-end reads to P-element, respectively.
mkdir p_alignment
echo "align "$read1" and $read2 to $refP"
echo "cmd: bowtie2 --local -q --no-unal -x $refP -U ${read1} -S p_alignment/read1_p.sam"
echo "cmd: bowtie2 --local -q --no-unal -x $refP -U ${read2} -S p_alignment/read2_p.sam"
bowtie2 --local -q --no-unal -x $refP -U $read1 -S p_alignment/read1_p.sam
bowtie2 --local -q --no-unal -x $refP -U $read2 -S p_alignment/read2_p.sam
echo "Aligning to reference genome is done\n"

# convert .sam to .bed
echo "convert .sam to .bed"
cd p_alignment
samtools view -bS read1_p.sam > read1_p.bam
bedtools bamtobed -i read1_p.bam -tag NM > read1_p.bed
samtools view -bS read2_p.sam > read2_p.bam
bedtools bamtobed -i read2_p.bam -tag NM > read2_p.bed
rm read1_p.bam read2_p.bam
cd ..
echo "conversion is done"
COMMENT

# pick reads mapped to P-elements and combine read1 and read2 (aligned to P-element) to a table
perl pick_p_mapped_reads.pl p_alignment/read1_p.bed p_alignment/read2_p.bed p_mapped_read.list
perl combine_paired_bed.pl p_alignment/read1_p.bed p_alignment/read2_p.bed p_mapped_read.list combined_paired_p.table
rm p_alignment/read1_p.bed p_alignment/read2_p.bed

# pick potential reads
perl pick_given_reads.pl -list p_mapped_read.list -i $read1 -out potential_read1.fq
perl pick_given_reads.pl -list p_mapped_read.list -i $read2 -out potential_read2.fq

# pick P-element-trimmed reads
perl pick_p_trimmed_seq.pl -mapped_list combined_paired_p.table -potential_read1 potential_read1.fq -potential_read2 potential_read2.fq -read1_sam p_alignment/read1_p.sam -read2_sam p_alignment/read2_p.sam -single_read1 single_read1.fq -single_read2 single_read2.fq -paired_read1 paired_read1.fq -paired_read2 paired_read2.fq

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
perl paired_breakpoint.pl combined_paired_p.table paired_dm6.bed counted_reads.txt paired_insertions.txt
perl single_breakpoint.pl combined_paired_p.table paired_insertions.txt single_dm6.bed counted_reads.txt paired_single_insertions.txt
perl combine_insertion_within_line.pl paired_single_insertions.txt combined_paired_single_insertions.txt
perl recruit_more_paired_reads.pl counted_reads.txt paired_dm6.bed combined_paired_single_insertions.txt combined_paired_p.table recruited_paired_insertions.txt
perl recruit_more_single_reads.pl counted_reads.txt recruited_paired_insertions.txt combined_paired_p.table single_dm6.bed $output

rm paired_dm6.sam single_dm6.sam paired_dm6.bed single_dm6.bed
rm paired_insertions.txt paired_single_insertions.txt combined_paired_single_insertions.txt recruited_paired_insertions.txt
rm counted_reads.txt combined_paired_p.table p_mapped_read.list
echo "...P-element insertion breakpoints are determined"


<<COMMENT
# estimate P-element insertion frequency
echo "...estimating insertion frequency begins..."

# build a reference containing constructed pseudo-genomes
echo "...build pseudogenomes..."
mkdir p_pseudo
perl make_P_pseudo_only.pl -f 500 -g $refGenome -p $refP -i $output -o p_pseudo/p_pseudo.fa
cd p_pseudo
cat p_pseudo.fa $refGenome > combined.fa
rm p_pseudo.fa
mv combined.fa p_pseudo.fa
bowtie2-build p_pseudo.fa p_pseudo.fa
cd ..
echo "...pseudo-genomes are built..."

# align paired-end reads to new reference
echo "...align to pseudo-genomes..."
bowtie2 -q --no-mixed --no-discordant --no-unal -x p_pseudo/p_pseudo.fa -1 $read1 -2 $read2 -S pseudo.sam
samtools view -bSf 0x2 -q 20 pseudo.sam > pseudo.bam
bedtools bamtobed -ed -i pseudo.bam > pseudo.bed
rm pseudo.bam
rm -r p_pseudo
echo "...aligning to pseudo-genomes is done..."

# calculate insertion frequency
perl freqCal_whole_genome_without_binSearch.pl -f 500 -i $output -b pseduo.bed -o temp_freq.txt
#rm $output
#mv temp_freq.txt $output
COMMENT


echo "...PIDFE is done..."
echo `date`



