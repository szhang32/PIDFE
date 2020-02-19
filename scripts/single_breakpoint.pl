#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
use Getopt::Long;

my %mapped_reads;
my $max_ed = 3;

#GetOptions;

my %table;
open IN, '<', $ARGV[0] or die "Cannot open the combined_paired_p.table file: $!";
while (my $line = <IN>) {
	chomp $line;
	my $read_id = (split /\t/, $line)[0];
	$table{$read_id} = $line;
}
close IN;

my %insert;
open IN, '<', $ARGV[1] or die "Cannot open the paired annnotated breakpoints: $!";
while (my $line = <IN>) {
    chomp $line;
    my ($key, $value) = split /\t/, $line;
    $insert{$key} = $value;
}
close IN;

open IN, '<', $ARGV[2] or die "Cannot open uniq_dm3.bed file: $!";
open OUT, '>>', $ARGV[3] or die "Cannot open the counted_reads.txt file: $!";
while (my $line = <IN>) {
    chomp $line;
	my ($chrom, $start, $end, $read, $edit_distance, $strand) = split /\t/, $line;
	my @bed = split /\t/, $table{$read};
	next unless $bed[1] eq "1p1";
	next if $edit_distance > $max_ed;

	my $direction = "";
    my $breakpoint = 0;
	if ($bed[3] - $bed[2] < $bed[7] - $bed[6]) { #read1 soft clipped
		next unless $bed[2] == 1 or $bed[3] == 2907; #no gap between P-element and adjacent genomic sequences
		if ($strand eq "+") { 
			$breakpoint = $end;
			if ($bed[5] eq "+") {
				$direction = "+";
			}
			elsif ($bed[5] eq "-") {
				$direction = "-";
			}
		}
		else {
			$breakpoint = $start+8;
			if ($bed[9] eq "+") {
                $direction = "+";
			}
			if ($bed[9] eq "-") {
                $direction = "-";
            }
		}
	}
	elsif ($bed[3] - $bed[2] > $bed[7] - $bed[6]) { #read2 soft clipped
		next unless $bed[6] == 1 or $bed[7] == 2907; #no gap between P-element and adjacent genomic sequences
		if ($strand eq "+") {
			$breakpoint = $end; 
			if ($bed[9] eq "+") {
				$direction = "+";
			}
			if ($bed[9] eq "-") {
				$direction = "-";
			}
		}
		else {
			$breakpoint = $start+8;
			if ($bed[5] eq "+") {
                $direction = "+";
            }
			if ($bed[5] eq "-") {
				$direction = "-";
			}
		}
	}	
	if ($direction ne ""){
        my $loc = $chrom.":".$breakpoint.":".$direction;
        $insert{$loc}++;
		print OUT "$read\n";
    }
}
close IN;
close OUT;

open OUT, '>', $ARGV[4] or die "Cannot open the output file: $!\n"; 
foreach my $key (sort keys %insert){
	print OUT  "$key\t$insert{$key}\n";
}
close OUT;

