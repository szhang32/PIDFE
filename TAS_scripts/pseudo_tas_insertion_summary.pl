#!/usr/bin/perl
use strict;
use warnings;

my ($strain, $infile, $outfile) = @ARGV;
open IN, '<', $infile or die;
open OUT, '>>', $outfile or die;
my $best_distance_insert = "";
my $best_distance_count = 0;
my $best_pseudo_insert = "NA";
my $best_pseudo_count = 0;
my $count = 0;
while (my $line = <IN>) {
	chomp $line;
	my ($insert, $distance_method, $pseudo_genome) = split /\t/, $line;
	next if $insert eq "ID";
	if ($distance_method > $best_distance_count) {
		$best_distance_insert = $insert;
		$best_distance_count = $distance_method 
	}	
	if ($pseudo_genome > $best_pseudo_count) {
		$best_pseudo_insert = $insert;
		$best_pseudo_count = $pseudo_genome
	}
	$count++;
}

if ($count == 0) {
	print OUT "$strain\t0\tNA\tNA\tNA\tNA\tNA\n";
}
elsif ($count == 1) {
	print OUT "$strain\t1\t$best_distance_insert\t$best_distance_count\t$best_pseudo_insert\t$best_pseudo_count\tNA\n";
}
elsif ($best_distance_insert eq $best_pseudo_insert) {
	print OUT "$strain\t$count\t$best_distance_insert\t$best_distance_count\t$best_pseudo_insert\t$best_pseudo_count\ttrue\n";
}
elsif ($best_distance_insert ne $best_pseudo_insert) {
	print OUT "$strain\t$count\t$best_distance_insert\t$best_distance_count\t$best_pseudo_insert\t$best_pseudo_count\tfalse\n";
}

close IN;
close OUT;
