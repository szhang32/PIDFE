# PIDFE
__PIDFE__, __*P*__-element **I**nsertion **D**etector and **F**requency **E**stimator, is a pipeline to detect *P*-element insertions and estimate their insertion frequencies from paired-end reads.

Copyright (c) 2019 Kelleher Lab at the University of Houston

Current version v1.0



## Installation
The following codes install PIDFE in home directory. Users can install it anywhere they want.

    cd ~
    git clone https://github.com/szhang32/PIDFE.git
    chmod 755 -R ~/PIDFE/scripts
    echo 'export PATH=$PATH:~/PIDFE/scripts'  >> ~/.bash_profile
    
## Software dependancies
PIDFE is designed to run on a high performance computering platform with Linux operating system.
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools](http://www.htslib.org/doc/samtools-1.2.html)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [BioPerl](https://bioperl.org)

## Usage
    sh PIDFE.sh [-i inDir] [-o output] <-g refGenome> <-p refP>
    
    Required arguments:
        <-g refGenome> defines the reference genome where *P*-elements reside.
        <-p refP> defines the consensus *P*-element sequence.
 
    Optional arguments:
        [-i inDir] is the directory containing paired-end read, i.e., *_R1.fastq and *_R2.fastq. Asterisk(*) can represent zero or any number of characters. Default: current working directoary  
        [-o output] is the output file name. Default: p_insertions.txt

## Contacts
If you find any bugs or have difficulties using PIDFE, please feel free to contact Shuo Zhang (shuozhang23@gmail.com).
