#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

###
#This Perl script is used to make a P element containing pseudo reference
#Usage: perl $0 [options]  std_ref insertion  pseudo
#Author: Shuo Zhang
#History: 2015-11-19
###

#read fasta from old reference to a hash
my %hash;
my %insertion;
my $in = Bio::SeqIO->new(-file => "$ARGV[0]" , '-format' => 'Fasta');
while( my $read = $in->next_seq()){
	my $name=$read->display_id;
	my $seq=$read->seq;
	$hash{$name}=$seq;
}

#read P element sequence from
my $p_seq; 
open IN,'<',$ARGV[1] or die "Cannot open the P element  file!\n";
while(my $line=<IN>){
	chomp $line;
	next if $line=~/>/;
	$p_seq .= $line;
}
close IN;
my $sense_seq=$p_seq;
$p_seq =~ tr /ATCG/TAGC/;
my $revComp=reverse($p_seq);

open INS,'<',$ARGV[2] or die "Cannot open the insertion file!\n";
open OUT,'>',$ARGV[3] or die "Cannot open the output file";
while (my $line=<INS>){
	chomp $line;
	my ($chrom, $pos, $strand) = split /:/, $line;
	print OUT ">$line\n";
	$pos = $pos - 8;
	my $pseudo="";
	my $begin=substr($hash{$chrom},$pos-500,500);
	my $end=substr($hash{$chrom},$pos+8,500);
	my $target_dup=substr($hash{$chrom},$pos,8);
	if ($strand eq "+"){
		$pseudo=$begin.$target_dup.$sense_seq.$target_dup.$end;
	}
	elsif ($strand eq "-"){
		$pseudo=$begin.$target_dup.$revComp.$target_dup.$end;
	}
	while (my $chunk = substr($pseudo, 0,50 ,"")) {
		print OUT  "$chunk\n";
	}
}
close INS;
close OUT;


