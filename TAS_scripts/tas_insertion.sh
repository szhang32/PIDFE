#!/bin/sh
# This shell script is used to detect P-element insertions in TAS (chr2R_TAS, chr3R_TAS, X_TAS and chr2L_TAS)
# 2016-11-08
# Shuo Zhang (shuozhang23@gmail.com)
# Erin Kelleher lab
# Program in Ecology and Evolution, Departement of Biology and Biochemistry, University of Houston

###
# cd to the directory containing all DGRP data: please run PIDF.sh first
# In the directory, create a directory named shell and put all scripts in the shell directory
# Please provide the path to non_tas_ref and tas_ref
###

cd /path/to/your_data/
### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

module add bowtie2/2.1.0
module add samtools/0.1.19
module add bedtools

non_tas_ref="/path/to/non_tas_ref"
tas_ref="/path/to/tas_ref"
for i in DGRP*
do 
	cd $i
	#align potential reads to non_tas and tas reference
	##align to non_tas first
	bowtie2 -q --un-conc paired_non_tas_unmapped.fq -x $non_tas_ref -1 paired_read1.fq -2 paired_read2.fq -S paired_non_tas.sam
	bowtie2 -q --un uniq_non_tas_unmapped.fq -x $non_tas_ref -U uniq_read1.fq,uniq_read2.fq -S uniq_non_tas.sam

	##align to tas reference
	bowtie2 -q --no-unal -a -x $tas_ref -1 paired_non_tas_unmapped.1.fq -2 paired_non_tas_unmapped.2.fq -S paired_tas.sam
	samtools view -XS paired_tas.sam > temp_paired_tas.sam
	bowtie2 -q  --no-unal -a -x $tas_ref -U uniq_non_tas_unmapped.fq -S uniq_tas.sam
	samtools view -XS uniq_tas.sam > formatted_uniq_tas.sam
	awk ' $2 ~ /U/' temp_paired_tas.sam >> formatted_uniq_tas.sam
	awk ' $2 ~ /P/' temp_paired_tas.sam > formatted_paired_tas.sam
	
	rm paired_non_tas_unmapped.?.fq
	rm uniq_non_tas_unmapped.fq

	#convert .sam to .bam
	samtools view -bS -o paired_tas.bam paired_tas.sam
	samtools view -bS -o uniq_tas.bam uniq_tas.sam 
	
	#rm paired_tas.sam uniq_tas.sam temp_paired_tas.sam
	perl ../shell/paired_distanceCal.pl  -i formatted_paired_tas.sam -o distance.txt
	perl ../shell/single_distanceCal.pl -i formatted_uniq_tas.sam -o distance.txt
	perl ../shell/distanceSummary.pl distance.txt $i >> ../distanceSummary.txt


    perl ../shell/reformat_distance.pl distance.txt reformatted_distance.txt
    #define the breakpoint
    perl ../shell/tas_paired_breakpoint.pl reformatted_distance.txt combined_paired_p.table formatted_paired_tas.sam tas_paired_insertion.txt
    perl ../shell/tas_single_breakpoint.pl reformatted_distance.txt combined_paired_p.table tas_paired_insertion.txt formatted_uniq_tas.sam tas_insertion.txt
    if [ -s combined_tas_insertion.txt ]
    then
        rm combined_tas_insertion.txt
    fi
    perl ../shell/combine_insertion_within_line.pl tas_insertion.txt combined_tas_insertion.txt
    perl ../shell/combine_XTAS_B_repeat.pl ../XTAS_repeat_coordinate_conversion.txt combined_tas_insertion.txt combined_XTAS_tas_insertion.txt
    rm tas_insertion.txt combined_tas_insertion.txt reformatted_distance.txt
    mv combined_XTAS_tas_insertion.txt combined_tas_insertion.txt

    #make psedugo
    rm -r tas_pseudo
    if [ -s tas_pseudo.sam ]
    then
        rm tas_pseudo.sam
    fi
    mkdir tas_pseudo
    cut -f 1 combined_tas_insertion.txt > temp.txt
    perl ../shell/pseudo_xtas.pl /project/kelleher/shuo/references/tas_ref/tas.fa /project/kelleher/shuo/references/p_element/p_element.fasta temp.txt tas_pseudo/tas_pseudo.fa
    rm temp.txt
    if [ -s tas_pseudo/tas_pseudo.fa ]
    then
        cd tas_pseudo
        bowtie2-build tas_pseudo.fa tas_pseudo.fa
    cd ..

    #pick potential read pairs
    grep -v 'ID' distance.txt  | cut -f 1 > temp.txt
    perl ../shell/pick_given_reads_v2.pl -list temp.txt -i potential_read1.fq -out tas_read1.fq
    perl ../shell/pick_given_reads_v2.pl -list temp.txt -i potential_read2.fq -out tas_read2.fq
    rm temp.txt

    #alignment
    bowtie2 -q --no-unal -x ./tas_pseudo/tas_pseudo.fa -1 tas_read1.fq -2 tas_read2.fq -S tas_pseudo.sam
    rm tas_read1.fq tas_read2.fq

    #samtools view -bSf 0x2 tas_pseudo.sam > temp.bam
    #bedtools bamtobed -i temp.bam > temp.bed
    #perl ../shell/count_tas_p.pl $i temp.bed

    samtools view -XSf 0x2 tas_pseudo.sam > formatted_tas_pseudo.sam
    rm tas_pseudo.sam

    rm pseudo_tas_insertion.txt
    perl ../shell/insertion_summary.pl combined_tas_insertion.txt formatted_tas_pseudo.sam  pseudo_tas_insertion.txt
    perl ../shell/coordinate_convertor.pl pseudo_tas_insertion.txt temp.txt
    rm pseudo_tas_insertion.txt formatted_tas_pseudo.sam
    mv temp.txt pseudo_tas_insertion.txt


    perl ../shell/pseudo_tas_insertion_summary.pl $i pseudo_tas_insertion.txt ../pseudo_tas_insertion.summary

    else
        echo "$i    NA    NA    NA    NA    NA    NA" >> ../pseudo_tas_insertion.summary
    fi
	cd ..
done
