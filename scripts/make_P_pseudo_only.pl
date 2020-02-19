#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

###
#This Perl script is used to make a P element containing pseudo reference
#Usage: perl $0 [options]  std_ref insertion  pseudo
#Author: Shuo Zhang
#History: 2019-06-21
###

my $flank_len = 500; # default flanking length is 500 bp
my $genome_ref;
my $p_ref;
my $insertions;
my $output;

GetOptions('f:i' => \$flank_len, 'g=s' => \$genome_ref, 'p=s' => \$p_ref, 'i=s' => \$insertions, 'o=s' => \$output);

#read fasta from old reference to a hash
my %hash;
my %insertion;
my $in = Bio::SeqIO->new(-file => "$genome_ref" , '-format' => 'Fasta');
while( my $read = $in->next_seq()){
	my $name=$read->display_id;
	my $seq=$read->seq;
	$hash{$name}=$seq;
}

#read P element sequence from
my $p_seq; 
open IN,'<',$p_ref or die "Cannot open the P element  file!\n";
while(my $line=<IN>){
	chomp $line;
	next if $line=~/>/;
	$p_seq .= $line;
}
close IN;
my $sense_seq=$p_seq;
$p_seq =~ tr /ATCG/TAGC/;
my $revComp=reverse($p_seq);

open INS,'<',$insertions or die "Cannot open the insertion file!\n";
open OUT,'>',$output or die "Cannot open the output file";
while (my $line=<INS>){
	chomp $line;
	next if $line=~/Chromosome/;
	my $ins = (split /\t/,$line)[0];
	my ($chrom, $pos, $strand) = split /:/,$ins;
	print OUT ">$ins:inserted\n";
	$pos = $pos - 8;
	my $pseudo="";
	my $begin=substr($hash{$chrom},$pos-$flank_len,$flank_len);
	my $end=substr($hash{$chrom},$pos+8,$flank_len);
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

	# get un_inserted contig
	#my $un_inserted = substr($hash{$chrom},$pos-$flank_len, $flank_len*2);
	#print OUT ">$ins:uninserted\n";
	#while (my $chunk = substr($un_inserted, 0,50 ,"")) {
		#print OUT  "$chunk\n";
	#}
}
close INS;
close OUT;




