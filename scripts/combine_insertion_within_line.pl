#!/usr/bin/perl
use strict;
use warnings;

my $step = 4;
open IN, '<', $ARGV[0] or die "Cannot open the input file: $!";
open OUT, '>', $ARGV[1] or die "Cannot open the output file: $!";
my $temp_strand = "";
my $max_count = 0;
my $total = 0;
my $temp_chrom = "";
my $temp_breakpoint = 0;
while (my $line = <IN>) {
	chomp $line;
	my ($insert, $count) = split /\t/, $line;
	my ($chrom, $breakpoint, $strand) = split /:/, $insert;	
	if ($temp_chrom ne $chrom or $temp_strand ne $strand) {
		if ($temp_chrom ne ""){
			my $loc = $temp_chrom.":".$temp_breakpoint.":".$temp_strand;
			print OUT "$loc\t$total\n";	
		}
		$temp_chrom = $chrom;
		$temp_breakpoint = $breakpoint;
		$temp_strand = $strand;
		$max_count = $count;
		$total = $count;
	}
	else {
		if (abs($breakpoint - $temp_breakpoint) > $step){
			my $loc = $temp_chrom.":".$temp_breakpoint.":".$temp_strand;
			print OUT "$loc\t$total\n";
			$temp_breakpoint = $breakpoint;	
			$temp_strand = $strand;
			$max_count = $count;
			$total = $count; 
		}
		else {
			if ($max_count < $count) {
				$max_count = $count;
				$temp_strand = $strand;
				$temp_breakpoint = $breakpoint;
			}
			$total += $count;
		}
	}
}
close IN;

my $loc = $temp_chrom.":".$temp_breakpoint.":".$temp_strand;
print OUT "$loc\t$total\n";

close OUT;

