# PIDFE
__PIDFE__, __*P*__-element **I**nsertion **D**etector and **F**requency **E**stimator, is a pipeline to detect *P*-element insertions and estimate their insertion frequencies from paired-end reads.

Current version v1.0

If you find any bugs or have difficulties using PIDFE, please feel free to contact Shuo Zhang (shuozhang23@gmail.com).


# Installation
The following codes install PIDFE in home directory. Users can install it anywhere they want.

    cd ~
    git clone https://github.com/szhang32/PIDFE.git
    echo 'export PATH=$PATH:~/PIDFE/scripts'  >> ~/.bash_profile
    
# Software dependencies
PIDFE is designed to run on a high performance computering platform with Linux operating system.
- [samtools](http://www.htslib.org/doc/samtools-1.2.html)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [Bio]
