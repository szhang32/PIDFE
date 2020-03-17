#!/usr/bin/perl -w 
use strict;

use Getopt::Long;

my $infile;
my $outfile;
my $cutoff = 3; #default cutoff
GetOptions ("in=s" => \$infile, "out=s" => \$outfile, "cutoff" => \$cutoff);
open IN,'<',$infile or die "Cannot open the sam file: $!";
open MAP,'>>',$outfile or die "Cannot open the output file: $!";
#print MAP "read_ID\tchr2R_TAS\tchr3R_TAS\tchrX_TAS\n";

my %hash;
my $line=""; #read id in question
while (my $align = <IN>){
	chomp $align;
	my @read = split /\t/,$align,12;
	my ($mismatch_read) = $read[11] =~ /XM:i:(\d+)/;
	my ($gap_read) = $read[11] =~ /XO:i:(\d+)/;
	#my $distance = $mismatch_read + $gap_read;
	my $distance = $read[11] =~ /NM:i:(\d+)/;
	next if $distance > $cutoff;
	if (($line ne $read[0])){	
		if (exists $hash{$line} and !((${$hash{$line}}[0]==100 and ${$hash{$line}}[1]==100 and ${$hash{$line}}[2]==100))) {
			print MAP "$line\t${$hash{$line}}[0]\t${$hash{$line}}[1]\t${$hash{$line}}[2]\n";
			delete $hash{$line};
		}
		$line = $read[0];
		#initialize
		${$hash{$line}}[0]=100; #chr2R_TAS
        	${$hash{$line}}[1]=100; #chr3R_TAS
        	${$hash{$line}}[2]=100; #chrX_TAS
	} 
	if (($read[2] eq "chrX_TAS") and ${$hash{$line}}[2] > $distance){
                ${$hash{$line}}[2]=$distance;
        }
        elsif (($read[2] eq "chr2R_TAS") and ${$hash{$line}}[0] > $distance){
                ${$hash{$line}}[0]=$distance;
        }
        elsif (($read[2] eq "chr3R_TAS") and ${$hash{$line}}[1] > $distance){
                ${$hash{$line}}[1]=$distance;
        }
}
#print the last read id
if ( (keys %hash)==1 and !(${$hash{$line}}[0]==100 and ${$hash{$line}}[1]==100 and ${$hash{$line}}[2]==100)){
	print MAP "$line\t${$hash{$line}}[0]\t${$hash{$line}}[1]\t${$hash{$line}}[2]\n";
}
close IN;
close MAP;

