#!/usr/bin/perl
use strict;
use warnings;

my ($alignment, $input, $output) = @ARGV;

my %align;
open IN, '<', $alignment or die;
while (my $line = <IN>) {
	chomp $line;
	my ($pos_abc, $pos_d) = split /\t/, $line;
	if ($pos_abc <= 2371) {
		next;
	}
	elsif ($pos_abc <= 4243) {
		next if $pos_d =~ /NA/;
		$align{$pos_d} = $pos_abc;
	}
	elsif ($pos_abc <= 6101) {
		if ($pos_d =~ /NA/) {
			$align{$pos_abc} = "NA";
		}
		else {
			$align{$pos_abc} = $align{$pos_d};
		}
	}
}
close IN;

my %breakpoints;
open my $in, '<', $ARGV[1] or die;
open my $out, '>', $ARGV[2] or die;
while (my $line = <$in>)  {
	chomp $line;
	if ($line =~ /X_TAS/) {	
		my ($loc, $reads) = split /\t/, $line;
		my ($pos, $strand) = (split /:/, $loc)[1,2];
		if (4244 <= $pos and $pos <= 7676) {
			if ($align{$pos} ne "NA") {
				my $temp = "chrX_TAS".":".$align{$pos}.":".$strand;
				$breakpoints{$temp} += $reads;
			}
			else {
				$breakpoints{$loc} = $reads
			}
		}
		else {
			$breakpoints{$loc} += $reads;
		}
	}	
}
close $in;


open IN, '<', $input or die;
while (my $line = <IN>)  {
    chomp $line;
    if ($line =~ /X_TAS/) {
		my ($loc, $reads) = split /\t/, $line;
        my ($pos, $strand) = (split /:/, $loc)[1,2];
		#if (exists $align{$pos} and $align{$pos} eq "NA") {
		#	print $out "$loc\t$reads\n";
		#	delete $breakpoints{$loc};
		#}
	}
	else {
		print $out "$line\n";
	}
}

foreach my $key (sort keys %breakpoints) {
	print $out "$key\t$breakpoints{$key}\n";
}
close IN;
close $out;
