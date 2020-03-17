#!/usr/bin/perl
use strict;
use warnings;

open IN, '<', $ARGV[0] or die;
open OUT, '>', $ARGV[1] or die;
while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /ID/){
		print OUT "$line\n";
		next;
	}
	my ($id, $distance, $pseudo) = split /\t/, $line;
	my ($chrom, $pos, $strand) = split /:/, $id;
	if ($chrom eq "chr2R_TAS") {
		$pos += 25258060 - 1;
	}
	elsif ($chrom eq "chr3R_TAS") {
		$pos += 32073015 - 3;
	}
	print OUT "$chrom:$pos:$strand\t$distance\t$pseudo\n";
}
close IN;
close OUT;
