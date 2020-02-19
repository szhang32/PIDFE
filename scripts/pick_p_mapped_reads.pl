#!/usr/bin/perl
use strict;
use warnings;

my ($read1, $read2, $output) = @ARGV;
my %res;

# P-element mapped reads in read1.fq
open IN, '<', $read1 or die;
while (my $line = <IN>) {
	chomp $line;
	my $id = (split /\t/, $line)[3];
	$res{$id} = 1;
}
close IN;

# P-element mapped reads in read2.fq
open IN, '<', $read2 or die;
while (my $line = <IN>) {
	chomp $line;
	my $id = (split /\t/, $line)[3];
	$res{$id} = 1;
}
close IN;

# print to output file
open OUT, '>', $output or die "Cannot open the output file: $!";
foreach my $key (sort keys %res) {
	print OUT "$key\n";
}
close OUT;
