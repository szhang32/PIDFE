#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
my %mapped_reads;
my $max_ed = 3;

my %table;
open IN, '<', $ARGV[0] or die "Cannot open the combined_paired_p.table file: $!";
while (my $line = <IN>) {
	chomp $line;
	my $read_id = (split /\t/, $line)[0];
	$table{$read_id} = $line;
}
close IN;


open IN, '<', $ARGV[1] or die "Cannot open bed file: $!";
open OUT, '>', $ARGV[2] or die "Cannot open the output file: $!";
my %insert;
while (defined (my $first  = <IN>) && defined (my $second = <IN>)) {
	chomp $first;
	chomp $second;
	my @read1=split /\t/, $first;
	my @read2=split /\t/, $second;
	next if $read1[4] > $max_ed or $read2[4] > $max_ed; #filter alignment with too many mismatches and gaps
	my $id1 = (split /\//, $read1[3])[0];
	my $id2 = (split /\//, $read2[3])[0];
	die unless $id1 eq $id2; #paired bed format requirement

	my $breakpoint = 0;
	my $direction = "";
	my @bed = split /\t/, $table{$id1};
	my ($chrom, $start, $end, $read, $mapq, $strand);
	if ($bed[1] eq "1p") { #read1 soft clipped
			if ($read1[3] =~ /\/1/) {
				($chrom, $start, $end, $read, $mapq, $strand) = @read1;
			}
			else {
				($chrom, $start, $end, $read, $mapq, $strand) = @read2;
			}
			next unless $bed[2] == 1 or $bed[3] == 2907; #no gap between P-element and adjacent genomic sequence
			if ($strand eq "-") { #read1 mapped to '-' strand of the reference genome
				$breakpoint = $end;
				if (($bed[3] - $bed[2] > 31) and $bed[5] eq "+") {  #the inverted repeat of P-element is 31bp
					$direction = "-";
				}
				elsif (($bed[3] - $bed[2] > 31) and $bed[5] eq "-") {
					$direction = "+";
				}
			}
			else { #read1 mapped to "+" strand of the reference genome
				$breakpoint = $start+8; # 8 bp target site duplication and bed format starts from 0
				if (($bed[3] - $bed[2] > 31) and $bed[5] eq "+") {
					$direction = "+";	
				}
				elsif (($bed[3] - $bed[2] > 31) and $bed[5] eq "-") {
					$direction= "-";
				}
			}
		}
	elsif ($bed[1] eq "p1") { #read2 soft clipped
			if ($read1[3] =~ /\/2/) {
				($chrom, $start, $end, $read, $mapq, $strand) = @read1;
			}
			else {
				($chrom, $start, $end, $read, $mapq, $strand) = @read2;
			}
			next unless $bed[6] == 1 or $bed[7] == 2907;  #no gap between P-element and adjacent genomic sequence
			if ($strand eq "-") { #read2 mapped to "-" strand of the reference genome
				$breakpoint = $end;
				if (($bed[7] - $bed[6] > 31) and $bed[9] eq "+") {  #the inverted repeat of P-element is 31bp
                    $direction = "-";
                }
                elsif (($bed[7] - $bed[6] > 31) and $bed[9] eq "-") {
                    $direction = "+";
                }
			}
			else { #read2 mapped to "+" strand of the reference genome
				$breakpoint = $start+8;
                if (($bed[7] - $bed[6] > 31) and $bed[9] eq "+") {
                    $direction = "+";
                }
                elsif (($bed[7] - $bed[6] > 31) and $bed[9] eq "-") {
                    $direction= "-";
                }
			}
	}
	if ($direction ne "") {
		my $loc = $chrom.":".$breakpoint.":".$direction;
		$insert{$loc}++;
		print OUT "$id1\n";
	}
}
close IN;
close OUT;

open OUT, '>', $ARGV[3] or die "Cannot open the output file: $!";
foreach my $key (sort keys %insert){
	print OUT "$key\t$insert{$key}\n";
}
close OUT;




