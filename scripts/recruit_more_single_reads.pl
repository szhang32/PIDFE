#!/usr/bin/perl
use strict;
use warnings;

my $step = 350;
my $pair_step = 3; #if two reads mapped to P-element, i.e. "1p1"
my $min_reads = 6;
my $max_ed = 3;


my %counted_reads;
open IN, '<', $ARGV[0] or die "Cannot open the counted_reads.txt file: $!";
while (my $line = <IN>) {
	chomp $line;
	$counted_reads{$line} = 1;
}
close IN;


my %insert;
open IN, '<', $ARGV[1] or die "Cannot open recruited_paired_insertions.txt file: $!\n";
while (my $line = <IN>) {
	chomp $line;
	my ($insertion, $count) = split /\t/, $line;
	$insert{$insertion} = $count;
}
close IN;

my %table;
open IN, '<', $ARGV[2] or die "Cannot open the combined_paired_p.table file: $!";
while (my $line = <IN>) {
	chomp $line;
	my $read_id = (split /\t/, $line)[0];
	$table{$read_id} = $line;
}
close IN;

open IN, '<', $ARGV[3] or die "Cannot open the mapping bed file: $!";
open OUT, '>>', $ARGV[0] or die "Cannot open the counted_reads.txt file: $!";
while (my $line = <IN> ) {
	chomp $line;
	my ($chrom, $start, $end, $read, $edit_distance, $strand) = split /\t/, $line;
	next if exists $counted_reads{$read};
	next if $edit_distance > $max_ed;

	my @bed = split /\t/, $table{$read};
	my $breakpoint = 0;
	my $direction = "";
	my $temp_site = "";
	my $distance = -1;
	
	if ($bed[1] eq "1p") {
		if ($strand eq "+") {
				$breakpoint = $end;
				if ($bed[5] eq "+") {
					$direction = "-";
				}
				elsif ($bed[5] eq "-") {
					$direction = "+";
				}
				#choose the closest insert breakpoint
				foreach my $key (keys %insert) {
					my ($insert_chrom, $insert_breakpoint, $insert_strand) = split /:/, $key;
					next unless $insert_chrom eq $chrom;
					next unless $insert_strand eq $direction;
					my $temp_distance = $insert_breakpoint - $breakpoint;
					if ($temp_distance > 0 and $temp_distance < $step) {
						if ($distance == -1 or $temp_distance < $distance) {
							$distance = $temp_distance;
							$temp_site = $key;
						}
					}
				}
				if ($distance != -1) {
					$insert{$temp_site}++;
					print OUT "$read\n";
					$counted_reads{$read} = 1;
				}
		}
		elsif ($strand eq "-") {
			$breakpoint = $start + 8;
			if ($bed[5] eq "+") {
				$direction = "+";
			}
			elsif ($bed[5] eq "-") {
				$direction = "-";
			}
			foreach my $key (keys %insert) {
				my ($insert_chrom, $insert_breakpoint, $insert_strand) = split /:/, $key;
				next unless $insert_chrom eq $chrom;
				next unless $insert_strand eq $direction;
				my $temp_distance = $breakpoint - $insert_breakpoint;
				if ($temp_distance > 0 and $temp_distance < $step) {
					if ($distance == -1 or $temp_distance < $distance) {
						$distance = $temp_distance;
						$temp_site = $key;
					}
				}
			}
			if ($distance != -1) {
				$insert{$temp_site}++;
				print OUT "$read\n";
				$counted_reads{$read} = 1;
			}
		}
	}
	elsif ($bed[1] eq "p1") {
		if ($strand eq "+") {
			$breakpoint = $end;
			if ($bed[9] eq "+") {
				$direction = "-";
			}
			elsif ($bed[9] eq "-") {
				$direction = "+";
			}
			foreach my $key (keys %insert) {
				my ($insert_chrom, $insert_breakpoint, $insert_strand) = split /:/, $key;
				next unless $insert_chrom eq $chrom;
				next unless $insert_strand eq $direction;
				my $temp_distance = $insert_breakpoint - $breakpoint;
				if ($temp_distance > 0 and $temp_distance < $step) {
					if ($distance == -1 or $temp_distance < $distance) {
						$distance = $temp_distance;
						$temp_site = $key;
					}
				}
			}
			if ($distance != -1) {
				$insert{$temp_site}++;
				print OUT "$read\n";
				$counted_reads{$read} = 1;
			}
		}
		elsif ($strand eq "-") {
			$breakpoint = $start + 8;
			if ($bed[9] eq "+") {
				$direction = "+";
			}
			elsif ($bed[9] eq "-") {
				$direction = "-";
			}
			foreach my $key (keys %insert) {
				my ($insert_chrom, $insert_breakpoint, $insert_strand) = split /:/, $key;
				next unless $insert_chrom eq $chrom;
				next unless $insert_strand eq $direction;
				my $temp_distance = $breakpoint - $insert_breakpoint;
				if ($temp_distance > 0 and $temp_distance < $step) {
					if ($distance == -1 or $temp_distance < $distance) {
						$distance = $temp_distance;
						$temp_site = $key;
					}
				}
			}
			if ($distance != -1) {
				$insert{$temp_site}++;
				print OUT "$read\n";
				$counted_reads{$read} = 1;
			}
		}			
	}
	elsif ($bed[1] eq "1p1") {
		if ($bed[3] - $bed[2] > $bed[7] - $bed[6]) { #read2 soft clipped
			if ($strand eq "+") {
				$breakpoint = $end;
		 		if ($bed[5] eq "+") {
					$direction = "-";
				}
				elsif ($bed[5] eq "-") {
					$direction = "+";
				}
			}
			if ($strand eq "-") {
				$breakpoint = $start + 8;
				if ($bed[5] eq "+") {
					$direction = "+";
				}
				elsif ($bed[5] eq "-") {
					$direction = "-";
				}
			}
		}
		else { #read1 soft clipped
			if ($strand eq "+") {
				$breakpoint = $end;
		 		if ($bed[9] eq "+") {
					$direction = "-";
				}
				elsif ($bed[9] eq "-") {
					$direction = "+";
				}
			}
			if ($strand eq "-") {
				$breakpoint = $start + 8;
				if ($bed[9] eq "+") {
					$direction = "+";
				}
				elsif ($bed[9] eq "-") {
					$direction = "-";
				}
			}
		}
		foreach my $key (sort keys %insert) {
			my ($insert_chrom, $insert_breakpoint, $insert_strand) = split /:/, $key;
			next unless $insert_chrom eq $chrom;
			next unless $insert_strand eq $direction;
			my $temp_distance = abs($insert_breakpoint - $breakpoint);
			if ($temp_distance <= $pair_step) {
				if ($distance == -1 or $temp_distance < $distance) {
					$distance = $temp_distance;
					$temp_site = $key;
				}
			}
		}
		if ($distance != -1) {
			$insert{$temp_site}++;
			print OUT "$read\n";
			$counted_reads{$read} = 1;
		}
	}
	else {
		die;
	}
}
close IN;
close OUT;

my $outfile = $ARGV[4];
open OUT, '>', $outfile or die "Cannot open the output file: $!";
foreach my $key (sort keys %insert){
	next if ($insert{$key} < $min_reads);
	print OUT "$key\t$insert{$key}\n";
	
}
close OUT;


