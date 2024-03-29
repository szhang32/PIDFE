# PIDFE
__PIDFE__, __*P*__-element **I**nsertion **D**etector and **F**requency **E**stimator, is a pipeline to detect *P*-element insertions and estimate their insertion frequencies from paired-end reads.

Copyright (c) 2019 Kelleher Lab at the University of Houston. If you used PIDFE in your study, please cite:
Zhang S, Pointer B, Kelleher E. 2020. Rapid evolution of piRNA-mediated silencing of an invading transposable element was driven by abundant de novo mutations. Genome Res. 30: 566-575

Current version v2.0

## How does PIDFE work
![](https://github.com/szhang32/PIDFE/blob/master/figures/PIDFE.png)
PIDFE identifies *P*-element insertions by detecting the split reads (read-1 in above figure), one part of which is mapped to *P*-element and the remaining part is mapped to the reference genome.

## Installation
The following codes install PIDFE in home directory. Users can install it anywhere they want by changing '~' to their desired directory.

    cd ~
    git clone https://github.com/szhang32/PIDFE.git
    chmod 755 -R ~/PIDFE/scripts
    echo 'export PATH=$PATH:~/PIDFE/scripts'  >> ~/.bash_profile
    
## Software dependancies
PIDFE is designed to run on a high performance computering platform with Linux operating system. The following softwares or packages are required to run PIDFE.
- [bowtie2/2.4.2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools/1.12](http://www.htslib.org/doc/samtools-1.2.html) 
- [bedtools/2.30.0](https://bedtools.readthedocs.io/en/latest/)
- [BioPerl/1.7.8](https://bioperl.org)

## Usage
    sh PIDFE.sh [-i inDir] [-o output] <-g refGenome> <-p refP>
    
    Required arguments:
        <-g refGenome> defines the reference genome where P-elements reside.
        <-p refP> defines the consensus P-element sequence.
 
    Optional arguments:
        [-i inDir] is the directory containing paired-end read, i.e., *_R1.fastq and *_R2.fastq. Asterisk(*) can represent zero or any number of characters. Default: current working directoary  
        [-o output] is the output file name. Default: p_insertions.txt

## Output format
The output file contains identified *P*-element insertion (1st column), the number of reads supporting each insertion (2nd column), and estimated frequencies (3rd column).
Here is an example:
| insertion | supporting_reads | frequency |
| ----------| -----------------| --------- |
| chr2L:12822166:+ | 27 | 0.844 |
| chr2L:20311155:- | 38 | 1 |
| chr2L:3632292:+ | 9	| 0.4285 |
| chr2L:3632292:- | 13 | 0.52 |

For each insertion, the inserted chromosome, insertion breakpoint, and orientation of the insertion are reported. For example, "chr2L:12822166:+" indicates that there is a sense *P*-element insertion at chromosome 2L at genomic coordinate 12822166. Sense ('+') represents the *P*-element insertion is in the same orientation as the plus reference genomic strand, whereas antisense ('-') represents the *P*-element insertion is in the same orientation as the minus reference genomic strand.

## Contacts
If you find any bugs or have difficulties using PIDFE, please feel free to contact Shuo Zhang (shuozhang23@gmail.com).
