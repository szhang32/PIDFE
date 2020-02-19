#!/usr/bin/perl
use strict;
use warnings;

###
#This script is used to recruit reads pairs that cannot determine P-element insertion orientation or there are mismatches
#around insertion breakpoint
###

die unless scalar @ARGV == 5;
my $max_ed = 3;
my $step = 3;


my %counted_reads;
open IN, '<', $ARGV[0] or die "Cannot open the counted_reads.txt file: $!";
while (my $line = <IN>) {
	chomp $line;
	$counted_reads{$line} = 1;
}
close IN;

my %insert;
open IN, '<', $ARGV[2] or die "Cannot open combined_paired_uniq_insertions.txt file: $!\n";
while (my $line = <IN>) {
	chomp $line;
	my ($insertion, $count) = split /\t/, $line;
	$insert{$insertion} = $count;
}
close IN;

my %table;
open IN, '<', $ARGV[3] or die "Cannot open the combined_paired_p.table file: $!";
while (my $line = <IN>) {
	chomp $line;
	my $read_id = (split /\t/, $line)[0];
	$table{$read_id} = $line;
}
close IN;


open IN, '<', $ARGV[1] or die "Cannot open the mapping bed file: $!";
open OUT, '>>', $ARGV[0] or die "Cannot open the counted_reads.txt file: $!";
while (defined (my $first  = <IN>) && defined (my $second = <IN>)) {
	chomp $first;
	chomp $second;
	my @read1=split /\t/, $first;
	my @read2=split /\t/, $second;
	next if $read1[4] > $max_ed or $read2[4] > $max_ed; #filter alignment with too many mismatches and gaps
	my $id1 = (split /\//, $read1[3])[0];
	my $id2 = (split /\//, $read2[3])[0];
	die unless $id1 eq $id2; #paired bed format requirement
	next if exists $counted_reads{$id1};
	
	my @bed = split /\t/, $table{$id1};	
	my $breakpoint = 0;
	my $direction = "";
	my ($chrom, $start, $end, $read, $edit_distance, $strand);
	if ($bed[1] eq "1p") { #read1 soft clipped
			if ($read1[3] =~ /\/1/) {
				($chrom, $start, $end, $read, $edit_distance, $strand) = @read1;
			}
			else {
				($chrom, $start, $end, $read, $edit_distance, $strand) = @read2;
			}
			#next unless $bed[2] == 1 or $bed[3] == 2907; #no gap between P-element and adjacent genomic sequence
			if ($strand eq "-") { #read1 mapped to '-' strand of the reference genome
				$breakpoint = $end;
				if ($bed[5] eq "+") {  #the inverted repeat of P-element is 31bp
					$direction = "-";
				}
				elsif ($bed[5] eq "-") {
					$direction = "+";
				}
			}
			else { #read1 mapped to "+" strand of the reference genome
				$breakpoint = $start+8; # 8 bp target site duplication and bed format starts from 0
				if ($bed[5] eq "+") {
					$direction = "+";	
				}
				elsif ($bed[5] eq "-") {
					$direction= "-";
				}
			}
		}
	elsif ($bed[1] eq "p1") { #read2 soft clipped
			if ($read1[3] =~ /\/2/) {
				($chrom, $start, $end, $read, $edit_distance, $strand) = @read1;
			}
			else {
				($chrom, $start, $end, $read, $edit_distance, $strand) = @read2;
			}
			#next unless $bed[6] == 1 or $bed[7] == 2907;  #no gap between P-element and adjacent genomic sequence
			if ($strand eq "-") { #read2 mapped to "-" strand of the reference genome
				$breakpoint = $end;
				if ($bed[9] eq "+") {  #the inverted repeat of P-element is 31bp
                    $direction = "-";
                }
                elsif ($bed[9] eq "-") {
                    $direction = "+";
                }
			}
			else { #read2 mapped to "+" strand of the reference genome
				$breakpoint = $start+8;
                if ($bed[9] eq "+") {
                    $direction = "+";
                }
                elsif ($bed[9] eq "-") {
                    $direction= "-";
                }
			}
	}

	my $distance = -1;
	my $temp_site = "";
	foreach my $key (sort keys %insert) {
		my ($insert_chrom, $insert_breakpoint, $insert_strand) = split /:/, $key;
		next unless $insert_chrom eq $chrom;
		#next unless $insert_strand eq $direction;
		my $temp_distance = abs($insert_breakpoint - $breakpoint);
		if ($temp_distance <= $step) {
			if ($distance == -1 or $temp_distance < $distance) {
				$distance = $temp_distance;
				$temp_site = $key;
			}
		}
	}
	if ($distance != -1) {
		$insert{$temp_site}++;
		print OUT "$id1\n";
		$counted_reads{$id1} = 1;
	}
}
close IN;
close OUT;

open OUT, '>', $ARGV[4] or die "Cannot open the output file: $!";
foreach my $key (sort keys %insert){
	print OUT "$key\t$insert{$key}\n";
}
close OUT;
print "\n";
