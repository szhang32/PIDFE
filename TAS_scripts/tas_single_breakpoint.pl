#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
my %mapped_reads;
open IN, '<', $ARGV[0] or die "Cannot open reformatted_distance.txt file: $!";
while (my $line = <IN>)
{
	chomp $line;		
	my ($read_id, $chrom, $distance) = split /\t/, $line;
	next if $read_id eq "read_ID";
	${$mapped_reads{$read_id}}[0] = $chrom;
	${$mapped_reads{$read_id}}[1] = $distance;
}
close IN;

my %table;
open IN, '<', $ARGV[1] or die "Cannot open the combined_paired_p.table file: $!";
while (my $line = <IN>) {
	chomp $line;
	my $read_id = (split /\t/, $line)[0];
	$table{$read_id} = $line;
}
close IN;

my %insertion;
open IN, '<', $ARGV[2] or die "Cannot open the tas_paired_insertion.txt file: $!";
while (my $line = <IN>) {
	chomp $line;
	my ($key, $value) = split /\t/, $line;
	$insertion{$key} = $value;
} 
close IN;


open IN, '<', $ARGV[3] or die "Cannot open formatted_uniq_tas.sam file: $!";
my ($first, $second);
my %temp;
my $count = 0;
my $temp_read_id = "";
while (my $line = <IN>) {
	chomp $line;
	my @sam = split /\t/, $line, 12;
	my $read_id = $sam[0];
	next unless (exists $mapped_reads{$read_id} and ${$mapped_reads{$read_id}}[0] eq $sam[2]);
	my ($mismatch) = $sam[11] =~ /XM:i:(\d+)/;
	my ($gap) = $sam[11] =~ /XO:i:(\d+)/;
	#my $distance = $mismatch + $gap;
	my $distance = $sam[11] =~ /NM:i:(\d+)/;
	next unless ${$mapped_reads{$read_id}}[1] == $distance;
	if ($temp_read_id ne $read_id) {
		if ($temp_read_id ne ""){
			foreach my $key (keys %temp) {
				if (exists $insertion{$key}) {
					$insertion{$key} += $temp{$key} / $count;
				}
				else {
					$insertion{$key} = $temp{$key} / $count;
				}
			}
		}	
		$temp_read_id = $read_id;
		$count = 0;
		%temp = ();
	}
	elsif ($temp_read_id eq $read_id or $temp_read_id eq "") {
		$temp_read_id = $read_id;
		my @bed = split /\t/, $table{$temp_read_id};
		my $breakpoint = 0;
		my $strand = "NA";
		my $aln_deletion = 0;
		if ($sam[1] =~ /p/ and $sam[1] =~ /1/ and $bed[1] eq "1p") {
			if ($sam[1] =~ /r/) {
				my (@cigar_m) = $sam[5] =~ /(\d+)M/g;
				if ($sam[5] =~ /D/) {
					my (@cigar_d) = $sam[5] =~ /(\d+)D/g;
					$aln_deletion = sum(@cigar_d);
				}
				my $aln_len = sum(@cigar_m);
				$breakpoint = $sam[3] + $aln_len + $aln_deletion - 1;
				if ($bed[3] - $bed[2] > 31 and $bed[5] eq "+") {
					$strand = "-";
				}
				elsif ($bed[3] - $bed[2] > 31 and $bed[5] eq "-") {
					$strand = "+";
				}
			}
			else {
				$breakpoint = $sam[3] + 7;
				if ($bed[3] - $bed[2] > 31 and $bed[5] eq "+") {
					$strand = "+";
				}
				elsif ($bed[3] - $bed[2] > 31 and $bed[5] eq "-") {
                	$strand = "-";
                }

			}
		}
		elsif ($sam[1] =~ /p/ and $sam[1] =~ /2/ and $bed[1] eq "p1") {
			if ($sam[1] =~ /r/) {
                my (@cigar_m) = $sam[5] =~ /(\d+)M/g;
                if ($sam[5] =~ /D/) {
                    my (@cigar_d) = $sam[5] =~ /(\d+)D/g;
                    $aln_deletion = sum(@cigar_d);
                }
                my $aln_len = sum(@cigar_m);
                $breakpoint = $sam[3] + $aln_len + $aln_deletion - 1;
                if ($bed[7] - $bed[6] > 31 and $bed[9] eq "+") {
                    $strand = "-";
                }
                elsif ($bed[7] - $bed[6] > 31 and $bed[9] eq "-") {
                    $strand = "+";
                }
            }
            else {
                $breakpoint = $sam[3] + 7;
                if ($bed[7] - $bed[6] > 31 and $bed[9] eq "+") {
                    $strand = "+";
                }
                elsif ($bed[7] - $bed[6] > 31 and $bed[9] eq "-") {
                    $strand = "-";
                }

            }
		}
		elsif ($bed[1] eq "1p1") {
			if ($bed[7] - $bed[6] > $bed[3] - $bed[2] and !($sam[1] =~ /r/)) {
				my (@cigar_m) = $sam[5] =~ /(\d+)M/g;
                if ($sam[5] =~ /D/) {
                    my (@cigar_d) = $sam[5] =~ /(\d+)D/g;
                    $aln_deletion = sum(@cigar_d);
                }
                my $aln_len = sum(@cigar_m);
                $breakpoint = $sam[3] + $aln_len + $aln_deletion - 1;
				if ($bed[9] eq "-") {
					$strand = "+";
				}
				elsif ($bed[9] eq "+") {
					$strand = "-"; 
				}
			}
			elsif ($bed[7] - $bed[6] < $bed[3] - $bed[2] and !($sam[1] =~ /r/)) {
				my (@cigar_m) = $sam[5] =~ /(\d+)M/g;
                if ($sam[5] =~ /D/) {
                    my (@cigar_d) = $sam[5] =~ /(\d+)D/g;
                    $aln_deletion = sum(@cigar_d);
                }
                my $aln_len = sum(@cigar_m);
                $breakpoint = $sam[3] + $aln_len + $aln_deletion - 1;
				if ($bed[5] eq "-") {
                    $strand = "+";
                }
                elsif ($bed[5] eq "+") {
                    $strand = "-";
                }
			}
			elsif ($bed[7] - $bed[6] < $bed[3] - $bed[2] and $sam[1] =~ /r/) {
				$breakpoint = $sam[3] + 7;
				if ($bed[5] eq "+") {
					$strand = "+";
				}
				elsif ($bed[5] eq "-") {
					$strand = "-";
				} 				
			}
			elsif ($bed[7] - $bed[6] > $bed[3] - $bed[2] and $sam[1] =~ /r/) {
				$breakpoint = $sam[3] + 7;
				if ($bed[9] eq "+") {
					$strand = "+";
				}
				elsif ($bed[9] eq "-") {            
                    $strand = "-";
                }
			}
			else {
				next;  #some reads $bed[7] - $bed[6] = $bed[3] - $bed[2] are ignored
			}
		}
		else {
			next;
		}
		my $loc = $sam[2].":".$breakpoint.":".$strand;
		$temp{$loc} = 1;
		$count++;
	}
}
close IN;

open OUT, '>', $ARGV[4] or die;
foreach my $key (sort keys %insertion) {
	print OUT "$key\t";
	printf OUT ("%.2f\n", $insertion{$key});
}


close OUT;
