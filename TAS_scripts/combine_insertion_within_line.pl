#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my @array = ("chr2R_TAS,+", "chr2R_TAS,-", "chr3R_TAS,+", "chr3R_TAS,-", "chrX_TAS,+", "chrX_TAS,-");
my $step = 3;
open IN, '<', $ARGV[0] or die "Cannot open the input file: $!";
while (my $line = <IN>) {
    chomp $line;
    my ($insert, $count) = split /\t/, $line;
    my ($chrom, $loc, $strand) = split /:/, $insert;
    #next unless $chrom eq "chrX_TAS";
	my $site = $chrom.":".$loc;
    ${$hash{$site}}{$strand} = $count;
}
close IN;


open OUT, '>', $ARGV[1] or die "Cannot open the output file: $!";
foreach my $element (@array){
my ($tas, $direction) = split /,/, $element;
open IN, '<', $ARGV[0] or die "Cannot open the input file: $!";
my $begin = 0;
my $potential = "";
my $max_count = 0;
my $total = 0;
my $flag = 0;
while (my $line = <IN>) {
	chomp $line;
	my ($insert, $count) = split /\t/, $line;
	my ($chrom, $loc, $strand) = split /:/, $insert;
	next unless ($chrom eq $tas and $strand eq $direction);
	my $site = $chrom.":".$loc;
	#combine "NA" to "-" or "+"
	if (exists ${$hash{$site}}{"NA"} and (!exists ${$hash{$site}}{"-"} or !exists ${$hash{$site}}{"+"})){
		$count += ${$hash{$site}}{"NA"};
	} 
	if (abs($loc - $begin) > $step) {
		if ($begin != 0) {
			print OUT "$potential\t$total\n";
		}
		$begin = $loc;
		$potential = $insert;
		$max_count = $count;
		$total = $count;
	}
	elsif (abs($loc - $begin) <= $step) {
		if ($count > $max_count){
			$potential = $insert;
			$max_count = $count;
		}
		$begin = $loc;
		$total += $count;
	}
}
close IN;
if ($begin != 0) {
	print OUT "$potential\t$total\n";
}
}

close OUT;

