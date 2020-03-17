#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my $mapq_cutoff = 20;
open IN, '<', $ARGV[0] or die;
while (my $line = <IN>) {
	chomp $line;
	my ($insert, $count) = split /\t/, $line;
	${$hash{$insert}}[0] = $count;
	${$hash{$insert}}[1] = 0;
}
close IN;

open IN, '<', $ARGV[1] or die;
my $first;
my $second;
while (defined ($first = <IN>) and defined ($second = <IN>)) {
	chomp $first;
	my @sam = split /\t/, $first;
	if ($sam[4] >= $mapq_cutoff and $sam[2] =~ /chr.*_TAS:(\d+):(\+|-)/) {
		${$hash{$sam[2]}}[1]++;
	}
}
close IN;

open OUT, '>', $ARGV[2] or die;
print OUT "ID\tdistance_method\tpseudo_genome\n";
foreach my $key (sort keys %hash){
	print OUT "$key\t${$hash{$key}}[0]\t${$hash{$key}}[1]\n";
}
close OUT;
