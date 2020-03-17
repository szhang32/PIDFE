#!/usr/bin/perl
use strict;
use warnings;

open IN, '<', $ARGV[0] or die "Cannot open distance.txt file: $!";
open OUT, '>', $ARGV[1] or die "Cannot open the reformatted.txt file: $!";
print OUT "read_ID\tchromosome\tdistance\n";
my $cutoff = 5; #maxmium mismatches and gaps
while (my $line = <IN>) {
	chomp $line;
	my @distance = split /\t/,$line;
	next if $distance[0] eq "read_ID";
	if ($distance[1] > $distance[3] and $distance[2] > $distance[3] and $distance[3] <= $cutoff) {
		print OUT "$distance[0]\tchrX_TAS\t$distance[3]\n";
	}
	elsif ($distance[1] < $distance[3] and $distance[1] < $distance[2] and $distance[1] <= $cutoff) {
		print OUT "$distance[0]\tchr2R_TAS\t$distance[1]\n";
	}
	elsif ($distance[2] < $distance[1] and $distance[2] < $distance[3] and $distance[2] <= $cutoff) {
		print OUT "$distance[0]\tchr3R_TAS\t$distance[2]\n";
	}
}
close IN;
close OUT;
