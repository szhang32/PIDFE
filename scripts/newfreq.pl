#!/usr/bin/perl -w 
use strict;
use Getopt::Long;

my $flank;
my $ins;
my $bed;
my $output;
my $bam;
GetOptions('flank=i' => \$flank, 'ins=s' => \$ins, 'bam=s' => \$bam, 'output=s' => \$output);

my $cross_start = $flank + 7;
my $cross_end = $flank + 2907 + 7 + 8;
my %hash;
my @sam;
my %all_ins;
my @sense;
my @anti;

open IN, '<', $ins or die "Cannot open the insertion file: $!";
(my $out = $ins) =~ s/insertion/freq/g;
open (OUTFILE, ">$out");
while (my $line=<IN>){
	chomp $line;
	#$line =~ s/>//g;
	my $insertion = (split /\t/, $line)[0];
	my $chrom = (split /:/, $insertion)[0];
	push(@{$all_ins{$chrom}}, $insertion);
	${$hash{$insertion}}[1]=0; ### reads supporting insertion
	#first reads forward strand
	@sense  = `samtools view -f 97 -q 20 sorted.new.pseudo.bam $insertion:inserted`;
	@anti  = `samtools view -f 161 -q 20 sorted.new.pseudo.bam $insertion:inserted`;
		
	for (@sense, @anti)
		{
		@sam = split /\t/,$_;
		if (($sam[3] < $cross_start and $sam[3] + $sam[8] > $cross_start) or ($sam[3] < $cross_end  and $sam[3] + $sam[8] > $cross_end)){ ${$hash{$insertion}}[1]++;}
		}
		
	(my $loc = $insertion) =~ s/:.$//g;
	(my $pos = $loc) =~ s/^.*?://g;
	${$hash{$insertion}}[0]= `samtools view -q 20 -c sorted.new.reference.bam $loc-$pos`;
	chomp(${$hash{$insertion}}[0]);
	my $freq = ${$hash{$insertion}}[1]/(${$hash{$insertion}}[0]+${$hash{$insertion}}[1]);
	###keep insertion only is more than 5 supporting reads
	if (${$hash{$insertion}}[1] > 5) {print OUTFILE $insertion."\t".${$hash{$insertion}}[0]."\t".${$hash{$insertion}}[1]."\t".$freq."\n";}
}
close IN;
close OUTFILE;

