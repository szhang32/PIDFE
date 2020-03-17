#!/usr/bin/perl -w 
use strict;

open IN,'<',$ARGV[0] or die "Cannot open the distance file: $!";
my $total_num=0;
my $chrX_TAS=0;
my $chr2R_TAS=0;
my $chr3R_TAS=0;
my $equal=0;
while (my $line=<IN>){
	chomp $line;
	next if $line =~ /read_ID/;
	next if $line eq "";
	$total_num++;
	my @distance = split /\t/,$line;
	if ($distance[1] > $distance[3] and $distance[2] > $distance[3]){
		$chrX_TAS++;
		#system("grep $distance[0] $file >> chrX_TAS_distance.sam ");
	}
	elsif ($distance[1] < $distance[3] and $distance[1] < $distance[2]) {
		$chr2R_TAS++;
	}
	elsif ($distance[2] < $distance[1] and $distance[2] < $distance[3]){
		$chr3R_TAS++;
	}
	else {
		$equal++;
	}
}
print "$ARGV[1]\t$total_num\t$chr2R_TAS\t$chr3R_TAS\t$chrX_TAS\t$equal\n";

