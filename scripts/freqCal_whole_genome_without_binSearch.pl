#!/usr/bin/perl -w 
use strict;
use Getopt::Long;

my $flank;
my $ins;
my $bed;
my $output;
GetOptions('flank=i' => \$flank, 'ins=s' => \$ins, 'bed=s' => \$bed, 'output=s' => \$output);

my $cross_start = $flank + 7;
my $cross_end = $flank + 2907 + 7 + 8;
my $max_ed = 3;
my %hash;
my %all_ins;

open IN, '<', $ins or die "Cannot open the insertion file: $!";
while (my $line=<IN>){
	chomp $line;
	#$line =~ s/>//g;
	my $insertion = (split /\t/, $line)[0];
	my $chrom = (split /:/, $insertion)[0];
	push(@{$all_ins{$chrom}}, $insertion);
	${$hash{$insertion}}[0]=0; #the number of reads that span the junction
	${$hash{$insertion}}[1]=0; #the number of reads that include P element
}
close IN;

open IN, '<', $bed or die "Cannot open the input bed file: $!";
while (defined (my $first = <IN>) && defined (my $second = <IN>)){
    chomp $first;
	chomp $second;
	my @read1=split /\t/,$first;
	my @read2=split /\t/,$second;
	next if $read1[4] > $max_ed or $read2[4] > $max_ed;
	if ($read1[0] =~ /inserted$/){
		my ($chrom, $pos, $strand, $inserted) = split /:/, $read1[0];
		my $insertion = $chrom.":".$pos.":".$strand;
		if ($read1[5] eq '+'){
			if (($read1[1] < $cross_start and $read2[2] > $cross_start) or ($read1[1] < $cross_end  and $read2[2] > $cross_end)){ ${$hash{$insertion}}[1]++;}
		}
		elsif ($read1[5] eq '-'){
			if (($read2[1] < $cross_start and $read1[2] > $cross_start) or ($read2[1] < $cross_end and $read1[2] > $cross_end)){ ${$hash{$insertion}}[1]++;}
		}
		else {
			die;
		}
	}
	else {
		next unless exists $all_ins{$read1[0]};
		my $start = 0;
		my $end = 0;
		if ($read1[5] eq '+') {
			$start = $read1[1];
			$end = $read2[2];
		}
		else {
			$start = $read2[1];
			$end = $read1[2]
		}
		foreach my $insertion (@{$all_ins{$read1[0]}}) {
			my ($chrom, $pos, $strand) = split /:/, $insertion;
			#my $flag = 0;
			if ($start < $pos and $pos < $end) {
				${$hash{$insertion}}[0]++;	
			}
			#if ($flag == 1) {
			#	if ($read1[5] eq '+') {
			#		last if $pos - $read2[2] > 2*$flank;
			#	}
			#	else {
			#		last if $pos - $read1[2] > 2*$flank;
			#	}
		}			
	}
}
close IN;

open OUT,'>',$output or die "Cannot open the output file: $!";
open IN,'<',$ins or die "Cannot open the input bed file: $!";
print OUT "insertion\tinserted_contig\tfraction\n";

while (my $line=<IN>){
	chomp $line;
	my $insertion = (split /\t/, $line)[0];
	if (${$hash{$insertion}}[0] + ${$hash{$insertion}}[1] > 0) {
		my $freq = ${$hash{$insertion}}[1] / (${$hash{$insertion}}[0] + ${$hash{$insertion}}[1]);
		print OUT "$insertion\t${$hash{$insertion}}[1]\t$freq\n";
	}
	else {
		print OUT "$insertion\t0\tNA\n";
	}
}
close OUT;
close IN;


