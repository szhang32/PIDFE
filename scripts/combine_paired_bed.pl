#!/usr/bin/perl 
use strict;
use warnings;

die unless @ARGV == 4;
my ($bed1, $bed2, $list, $outfile) = @ARGV;

my %read1_bed;
open BED, '<', $bed1 or die "Cannot open the bed1 file: $!";
while (my $line = <BED>){
	chomp $line;
	my $id = (split /\t/, $line)[3];
	$read1_bed{$id} = $line;
}
close BED;

my %read2_bed;
open BED, '<', $bed2 or die "Cannot open the bed2 file: $!";
while (my $line = <BED>){
    chomp $line;
    my $id = (split /\t/, $line)[3];
    $read2_bed{$id} = $line;
}
close BED;


open IN, '<', $list or die "Cannot open the mapped list file: $!";
open OUT, '>', $outfile or die "Cannot open the output file: $!";
print OUT "read_id\tcategory\tread1_start\tread1_end\tread1_edit_distance\tread1_strand\tread2_start\tread2_end\tread2_edit_distance\tread2_strand\n";
while (my $id = <IN>){
	chomp $id;
	if (exists $read1_bed{$id} and ! (exists $read2_bed{$id})){
		my @read1 = split /\t/, $read1_bed{$id};	
		my $read1_start = $read1[1] + 1;
		print OUT "$id\t1p\t$read1_start\t$read1[2]\t$read1[4]\t$read1[5]\t\*\t\*\t\*\t\*\n"; #1p means only read1 aligned to P-element
		delete $read1_bed{$id};
	}
	elsif (exists $read2_bed{$id} and ! (exists $read1_bed{$id})){
		my @read2 = split /\t/, $read2_bed{$id};
		my $read2_start = $read2[1] + 1;
		print OUT "$id\tp1\t\*\t\*\t\*\t\*\t$read2_start\t$read2[2]\t$read2[4]\t$read2[5]\n"; #p1 means only read2 aligned to P-element
		delete $read2_bed{$id};
	}
	else {
		my @read1 = split /\t/, $read1_bed{$id};
		my @read2 = split /\t/, $read2_bed{$id};
		if ($read1[5] eq $read2[5]) { #read1 and read2 mapped to same strand of P-element, one of them should be corrected due to 5' and 3' inverted repeat
			if ($read1[2] - $read1[1] > $read2[2] - $read2[1]) {
				#correct strand
				if ($read1[5] eq "+") {
					$read2[5] = "-";
				}
				else {
					$read2[5] = "+";
				}

				if ($read1[1] > 1000) { #read2 should be mapped to 3' end
					#die unless $read2[1] < 3;
					my $temp = $read2[1];
					$read2[1] = 2907 - ($read2[2] - $read2[1]);
					$read2[2] = 2907 - $temp;
				}
				else { #read2 should be mapped to 5' end
					#die unless (2907 - $read2[2]) < 3;
					my $temp = $read2[2];
					$read2[2] = $read2[2] - $read2[1];
					$read2[1] = 2907 - $temp;
				}	
			}
			else {
				#correct strand
				if ($read2[5] eq "+") {
					$read1[5] = "-";
				}
				else {
					$read1[5] = "+";
				}

				if ($read2[1] > 1000) { #read1 should be mapped to 3' end
					my $temp = $read1[1];
					$read1[1] = 2907 - ($read1[2] - $read1[1]);
					$read1[2] =	2907 - $temp;
				}
				else { #read1 should be mapped to 5' end
					my $temp = $read1[2];
					$read1[2] = $read1[2] - $read1[1];
					$read1[1] = 2907 - $temp;
				}
			}
		}
		my $read1_start = $read1[1] + 1;
		my $read2_start = $read2[1] + 1;
		print OUT "$id\t1p1\t$read1_start\t$read1[2]\t$read1[4]\t$read1[5]\t$read2_start\t$read2[2]\t$read2[4]\t$read2[5]\n"; #1p1 means both read1 and read2 aligned to P-element
		delete $read1_bed{$id}; 
		delete $read2_bed{$id};
	}
}
close IN;
die unless keys %read1_bed == 0 and keys %read2_bed == 0;

