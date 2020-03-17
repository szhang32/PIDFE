#!/usr/bin/perl -w 
use strict;
use Getopt::Long;

my $infile;
my $outfile;
my $cutoff = 3; #default cutoff
GetOptions ("in=s" => \$infile, "out=s" => \$outfile, "cutoff" => \$cutoff);
open IN,'<',$infile or die "Cannot open the sam file: $!";
open MAP,'>',$outfile or die "Cannot open the output file: $!";
print MAP "read_ID\tchr2R_TAS\tchr3R_TAS\tchrX_TAS\n";
my %hash;
my ($first, $second);
my $line=""; #read id in question
while (defined ($first = <IN>) && defined ($second = <IN>)){
	chomp $first;
	chomp $second;
	my @read1=split /\t/,$first,12;
	my @read2=split /\t/,$second,12;
	my ($mismatch_read1)=$read1[11]=~/XM:i:(\d+)/;
	my ($mismatch_read2)=$read2[11]=~/XM:i:(\d+)/;
	my ($gap_read1)=$read1[11]=~/XO:i:(\d+)/;
	my ($gap_read2)=$read2[11]=~/XO:i:(\d+)/;
	#my $distance=$mismatch_read1+$mismatch_read2+$gap_read1+$gap_read2;
	#next if $distance > 5;
	my ($distance1) = $read1[11] =~ /NM:i:(\d+)/;
	my ($distance2) = $read2[11] =~ /NM:i:(\d+)/;
	next if $distance1 > $cutoff or $distance2 > $cutoff;
	my $distance = $distance1 + $distance2;
	if (($line ne $read1[0])){
		if (exists $hash{$line} and !((${$hash{$line}}[0]==100 and ${$hash{$line}}[1]==100 and ${$hash{$line}}[2]==100))) {
			print MAP "$line\t${$hash{$line}}[0]\t${$hash{$line}}[1]\t${$hash{$line}}[2]\n";
			delete $hash{$line};
		}
		$line = $read1[0];
		#initialize
		${$hash{$line}}[0]=100; #chr2R_TAS
        	${$hash{$line}}[1]=100; #chr3R_TAS
        	${$hash{$line}}[2]=100; #chrX_TAS
	} 
	if (($read1[2] eq "chrX_TAS") and ${$hash{$line}}[2] > $distance){
                ${$hash{$line}}[2]=$distance;
        }
        elsif (($read1[2] eq "chr2R_TAS") and ${$hash{$line}}[0] > $distance){
                ${$hash{$line}}[0]=$distance;
        }
        elsif (($read1[2] eq "chr3R_TAS") and ${$hash{$line}}[1] > $distance){
                ${$hash{$line}}[1]=$distance;
        }
}
#print the last read id
if (keys %hash == 1 and !(${$hash{$line}}[0]==100 and ${$hash{$line}}[1]==100 and ${$hash{$line}}[2]==100)){
	print MAP "$line\t${$hash{$line}}[0]\t${$hash{$line}}[1]\t${$hash{$line}}[2]\n";
}
close IN;
close MAP;

