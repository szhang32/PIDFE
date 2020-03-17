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

open IN, '<', $ARGV[2] or die "Cannot open formatted_paired_tas.sam file: $!";
my ($first, $second);
my %insertion;
my %temp;
my $count = 0;
my $temp_read_id = "";
while (defined ($first = <IN>) && defined ($second = <IN>)) {
	chomp $first;
	chomp $second;
	my @read1=split /\t/,$first,12;
	my $read_id = $read1[0];
	next unless (exists $mapped_reads{$read_id} and ${$mapped_reads{$read_id}}[0] eq $read1[2]);
	my @read2=split /\t/,$second,12;
	my ($mismatch_read1)=$read1[11]=~/XM:i:(\d+)/;
	my ($mismatch_read2)=$read2[11]=~/XM:i:(\d+)/;
	my ($gap_read1)=$read1[11]=~/XO:i:(\d+)/;
	my ($gap_read2)=$read2[11]=~/XO:i:(\d+)/;
	#my $distance=$mismatch_read1+$mismatch_read2+$gap_read1+$gap_read2;
	my ($distance1) = $read1[11] =~ /NM:i:(\d+)/;
	my ($distance2) = $read2[11] =~ /NM:i:(\d+)/;
	my $distance = $distance1 + $distance2;
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
		if ($bed[1] eq "1p") { #read1 soft clipped
			if ($read1[1] =~ /r/) {  #read1 mapped to '-' strand of the reference genome
				my (@cigar_m) = $read1[5] =~ /(\d+)M/g;
				if ($read1[5] =~ /D/) {
					my (@cigar_d) = $read1[5] =~ /(\d+)D/g;
					$aln_deletion = sum(@cigar_d);
				}
				my $aln_len = sum(@cigar_m);
				$breakpoint = $read1[3] + $aln_len + $aln_deletion - 1;
				if ($bed[3] - $bed[2] > 31 and $bed[5] eq "+") {
					$strand = "-";
				}
				elsif ($bed[3] - $bed[2] > 31 and $bed[5] eq "-") {
					$strand = "+";
				}
			}
			else {  #read1 mapped to "+" strand of the reference genome
				$breakpoint = $read1[3] + 7;
				if ($bed[3] - $bed[2] > 31 and $bed[5] eq "+") {
					$strand = "+";
				}
				elsif ($bed[3] - $bed[2] > 31 and $bed[5] eq "-") {
                	$strand = "-";
                }
			}
		}
		elsif ($bed[1] eq "p1") { #read2 soft clipped
			if ($read2[1] =~ /r/) { #read2 mapped to "-" strand of the reference genome
				my (@cigar_m) = $read2[5] =~ /(\d+)M/g;
				if ($read1[5] =~ /D/) {
                    my (@cigar_d) = $read1[5] =~ /(\d+)D/g;
                    $aln_deletion = sum(@cigar_d);
                }
				my $aln_len = sum(@cigar_m);
                $breakpoint = $read2[3] + $aln_len + $aln_deletion - 1;
				if ($bed[7] - $bed[6] > 31 and $bed[9] eq "+") { 
					$strand = "-";
				}
				elsif ($bed[7] - $bed[6] > 31 and $bed[9] eq "-") {
					$strand = "+";
				}
			}
			else {
				$breakpoint = $read2[3] + 7;
				if ($bed[7] - $bed[6] > 31 and $bed[9] eq "+") { 
                    $strand = "+";
                }
                elsif ($bed[7] - $bed[6] > 31 and $bed[9] eq "-") {
                    $strand = "-";
                }
			}
		}
		else {
			next;
		}
		my $loc = $read1[2].":".$breakpoint.":".$strand;
		$temp{$loc} = 1;	
		$count++;	
	}
}
close IN;

open OUT, '>', $ARGV[3] or die "Cannot open the output file: $!";
foreach my $key (keys %insertion) {
	print OUT "$key\t";
	printf OUT ("%.2f\n", $insertion{$key});
}
close OUT;

