#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

###
##This perl script is usage to pick up read pair according to a list
##Usage: perl $0 [options] -read_list mapped.read.list input.fq out.fq
##Author: Shuo Zhang
##History: 2017-1-25
####
my $opposite = ''; #default print out the read in the list
my $list; 
my $infile;
my $outfile;
GetOptions('opposite!' => \$opposite, 'list=s' => \$list, 'infile=s' => \$infile, 'outfile=s' => \$outfile) or die ("Error in command line arguments\n");


#read p mapped read to a hash 
open LIST, '<', $list or die "Cannot open the p_mapped reads list: $!";
my %read_list;
while (my $line = <LIST>){
	chomp $line;
	$line = "\@".$line;
	$read_list{$line} = 1; #insert the read id to read_list
}
close LIST;

open IN, '<', $infile or die "Cannot open the input fastq file: $!";
open OUT, '>', $outfile or die "Cannot open the output fastq file: $!";
my ($read_id, $seq, $read_id_opt, $quality);
while (defined ($read_id = <IN>) and defined ($seq = <IN>) and defined ($read_id_opt = <IN>) and defined ($quality = <IN>)){
	chomp $read_id;
	chomp $seq;
	chomp $read_id_opt;
	chomp $quality;
	die unless ($read_id =~ /^\@/ and $read_id_opt =~ /^\+/); #read_id should begin with '@' and optional id should begin with '+'. This line is a checkpoint
	my $primary_id = (split / /, $read_id)[0];
	if (exists $read_list{$primary_id} and !$opposite){	
		print OUT "$primary_id\n";
		print OUT "$seq\n";
		print OUT "+\n";
		print OUT "$quality\n"; 
	} 
	elsif (!(exists $read_list{$primary_id}) and $opposite){
		print OUT "$primary_id\n";
        print OUT "$seq\n";
        print OUT "+\n";
        print OUT "$quality\n";
	}
}
close IN;
close OUT;
