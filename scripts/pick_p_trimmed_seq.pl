#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
use Getopt::Long;

###
#This perl script is used to extract reads mapped to P-element (one mismatch and one gap open are allowed)
#Usage: perl $0 out.sam read1.fq r1.fq read2.fq r2.fq
#History: 2015-09-08
###
my $min_len=30;
my $mapped_list;
my $read1_sam;
my $read2_sam;
my $potential_read1;
my $potential_read2;
my $single_read1;
my $single_read2;
my $paired_read1;
my $paired_read2;
GetOptions ('potential_read1=s' => \$potential_read1, 'potential_read2=s' => \$potential_read2, 'min:i'=>\$min_len, 'mapped_list=s' => \$mapped_list, 'read1_sam=s' => \$read1_sam,'read2_sam=s' => \$read2_sam, 'single_read1=s' => \$single_read1, 'single_read2=s' => \$single_read2, 'paired_read1=s' => \$paired_read1, 'paired_read2=s' => \$paired_read2) or die ("Error in command line arguments\n");

#sub function to obtain reverse complementary DNA sequence
sub reverse_complement {
	my $seq = shift;
	my $revcomp=reverse($seq);
	$revcomp=~tr/ATCG/TAGC/;
	return $revcomp;
}

#store read1 that is aligned to P-element to a hash %align1
my %align1;
open SAM1,'<',$read1_sam or die "Cannot open the read1_sam file: $!";
while (my $line1=<SAM1>){
	chomp $line1;
	next if $line1 =~ /^\@/;
	my @sam1=split /\t/,$line1,2;
	$align1{$sam1[0]}=$sam1[1];
}
close SAM1;

#store read2 that is aligned to P-element to a hash %align2
my %align2;
open SAM2,'<',$read2_sam or die "Cannot open the read2_sam file: $!";
while (my $line2=<SAM2>){
	chomp $line2;
	next if $line2 =~ /^\@/;
	my @sam2=split /\t/,$line2,2;
	$align2{$sam2[0]}=$sam2[1];
}
close SAM2;

#store potential paired reads to a hash %po_r1;
my ($read_id, $seq, $read_id_opt, $quality);
my %po_r1;
open IN,'<', $potential_read1 or die "Cannot open the potentail_read1.fq file: $!";
while (defined ($read_id = <IN>) and defined ($seq = <IN>) and defined ($read_id_opt = <IN>) and defined ($quality = <IN>)){
	chomp $read_id;
	chomp $seq;
	chomp $read_id_opt;
	chomp $quality;
	die unless ($read_id =~ /^\@/ and $read_id_opt =~ /^\+/); #read_id should begin with '@' and optional id should begin with '+'. This line is a checkpoint
	my $value = $seq."\t".$read_id_opt."\t".$quality;
	$po_r1{$read_id} = $value;
} 
close IN;

my %po_r2;
open IN,'<', $potential_read2 or die "Cannot open the potentail_read2.fq file: $!";
while (defined ($read_id = <IN>) and defined ($seq = <IN>) and defined ($read_id_opt = <IN>) and defined ($quality = <IN>)){
    chomp $read_id;
    chomp $seq;
    chomp $read_id_opt;
    chomp $quality;
    die unless ($read_id =~ /^\@/ and $read_id_opt =~ /^\+/); #read_id should begin with '@' and optional id should begin with '+'. This line is a checkpoint
	my $value = $seq."\t".$read_id_opt."\t".$quality;
    $po_r2{$read_id} = $value;
}
close IN;


open U_R1,'>',$single_read1 or die "Cannot open the single_read1.fq file: $!\n";
open U_R2,'>',$single_read2 or die "Cannot open the single_read2.fq file: $!\n";
open PAIR_R1,'>',$paired_read1 or die "Cannot open the paired_read1.fq file: $!\n";
open PAIR_R2,'>',$paired_read2 or die "Cannot open the paired_read2.fq file: $!\n";
open IN,'<',$mapped_list or die "Cannot open the mapped_list: $!";
while (my $line=<IN>){
	chomp $line;
	my @bed = split /\t/, $line;
	next if $bed[0] eq "read_id";
	my $id="\@$bed[0]";
	if ($bed[1] eq "1p"){ # read1 mapped to p element, read2 not
		my @sam1=split /\t/, $align1{$bed[0]};
		#die unless exists $po_r2{$id};
		if (! exists $po_r2{$id}) {
			print "$id is not existing!\n";
			next;
		}
		my $whole_read=$po_r2{$id};
		my @array=split /\t/, $whole_read;
		if ($sam1[4]=~/^(\d+)S(\d+)M/ and $1>=$min_len){ #left soft-clipped
			my $seq = substr $sam1[8],0,$1;
			my $quality=substr $sam1[9],0,$1;
            if ($sam1[0] & 16){ # the read is mapped to reverse strand
				my $revcomp=reverse_complement($seq);
				my $revcomp_quality=reverse($quality);
				print PAIR_R1 "$id\n$revcomp\n\+\n$revcomp_quality\n";
			}
			else {
				print PAIR_R1 "$id\n$seq\n\+\n$quality\n";
			}
			print PAIR_R2 "$id\n$array[0]\n$array[1]\n$array[2]\n";
		}		
		elsif ($sam1[4]=~/(\d+)M(\d+)S$/ and $2 >=$min_len){ #right soft-clipped
			my $seq=substr $sam1[8],-$2;
			my $quality=substr $sam1[9],-$2;
			if ($sam1[0] & 16){
				my $revcomp=reverse_complement($seq);
				my $revcomp_quality=reverse($quality);
				print PAIR_R1 "$id\n$revcomp\n\+\n$revcomp_quality\n";
			}
			else {
				print PAIR_R1 "$id\n$seq\n\+\n$quality\n";
			}
			print PAIR_R2 "$id\n$array[0]\n$array[1]\n$array[2]\n";
		}
		else {
			print U_R2 "$id\n$array[0]\n$array[1]\n$array[2]\n";
		}
		delete $align1{$bed[0]};
	}
	elsif ($bed[1] eq "p1"){ #read2 mapped to p element, read1 not
		my @sam2=split /\t/, $align2{$bed[0]};
		#die unless exists $po_r1{$id};
        if (! exists $po_r1{$id}) {
			print "$id is not existing!\n";
			next;
		}
		my $whole_read=$po_r1{$id};
        my @array=split /\t/, $whole_read;
        if ($sam2[4]=~/^(\d+)S(\d+)M/ and $1>=$min_len){ #left soft-clipped
            my $seq = substr $sam2[8],0,$1;
            my $quality=substr $sam2[9],0,$1;
            if ($sam2[0] & 16){
                my $revcomp=reverse_complement($seq);
                my $revcomp_quality=reverse($quality);
                print PAIR_R2 "$id\n$revcomp\n\+\n$revcomp_quality\n";
            }
            else {
                print PAIR_R2 "$id\n$seq\n\+\n$quality\n";
            }
            print PAIR_R1 "$id\n$array[0]\n$array[1]\n$array[2]\n";
        }
        elsif ($sam2[4]=~/(\d+)M(\d+)S$/ and $2 >=$min_len){ #right soft-clipped
            my $seq=substr $sam2[8],-$2;
            my $quality=substr $sam2[9],-$2;
            if ($sam2[0] & 16){
                my $revcomp=reverse_complement($seq);
                my $revcomp_quality=reverse($quality);
                print PAIR_R2 "$id\n$revcomp\n\+\n$revcomp_quality\n";
            }
            else {
                print PAIR_R2 "$id\n$seq\n\+\n$quality\n";
            }
            print PAIR_R1 "$id\n$array[0]\n$array[1]\n$array[2]\n";
        }
        else {
            print U_R1 "$id\n$array[0]\n$array[1]\n$array[2]\n";
        }
		delete $align2{$bed[0]};
	}
	elsif ($bed[1] eq "1p1"){ ##both reads mapped to p element
		my @sam1=split /\t/,$align1{$bed[0]};
		my $match_len1 = $bed[3] - $bed[2] + 1;
		my @sam2=split /\t/,$align2{$bed[0]};
		my $match_len2 = $bed[7] - $bed[6] + 1;
		if ($match_len1 <= $match_len2) { #read1 soft-clipped
			if ($sam1[4] =~ /(\d+)S(\d+)M/ and $1 >= $min_len){ #read1 left soft-clippped
				my $seq = substr $sam1[8],0,$1;
				my $quality=substr $sam1[9],0,$1;
				if ($sam1[0] & 16){
					my $revcomp=reverse_complement($seq);
					my $revcomp_quality=reverse($quality);
					print U_R1 "$id\n$revcomp\n\+\n$revcomp_quality\n";
				}
				else {
					print U_R1  "$id\n$seq\n\+\n$quality\n";
				}
			}
			elsif ($sam1[4] =~ /(\d+)M(\d+)S$/ and $2 >= $min_len) { #read1 right soft-clipped
				my $seq=substr $sam1[8],-$2;
				my $quality=substr $sam1[9],-$2;
				if ($sam1[0] & 16) {
					my $revcomp=reverse_complement($seq);
					my $revcomp_quality=reverse($quality);
					print U_R1 "$id\n$revcomp\n\+\n$revcomp_quality\n";
				}
				else {
					print U_R1 "$id\n$seq\n\+\n$quality\n";
				}
			}
		}
		else { #read2 soft-clipped
			if ($sam2[4] =~ /(\d+)S(\d+)M/ and $1 >= $min_len) {  #read2 left soft-clipped
				my $seq = substr $sam2[8],0,$1;
				my $quality=substr $sam2[9],0,$1;
				if ($sam2[0] & 16){
					my $revcomp=reverse_complement($seq);
                    my $revcomp_quality=reverse($quality);
					print U_R2 "$id\n$revcomp\n\+\n$revcomp_quality\n";
				}
				else {
					print U_R2  "$id\n$seq\n\+\n$quality\n";
				}
			}
			elsif ($sam2[4] =~ /(\d+)M(\d+)S$/ and $2 >= $min_len) { #read2 right soft-clipped
				my $seq=substr $sam2[8],-$2;
                my $quality=substr $sam2[9],-$2;
				if ($sam2[0] & 16){
	                my $revcomp=reverse_complement($seq);
                    my $revcomp_quality=reverse($quality);
					print U_R2 "$id\n$revcomp\n\+\n$revcomp_quality\n";
				}	
				else {
					print U_R2  "$id\n$seq\n\+\n$quality\n";
				}
			}
		}
		delete $align1{$bed[0]};
		delete $align2{$bed[0]};
	}
	else {
		die;
	}
	delete $po_r1{$id};
	delete $po_r2{$id};
}
close IN;
close U_R1;
close U_R2;
close PAIR_R1;
close PAIR_R2;
die unless keys %align1 == 0 and keys %align2 == 0 and keys %po_r1 == 0 and keys %po_r2 == 0;
print "The end!\n";

