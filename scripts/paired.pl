#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
my $file;
GetOptions ( "file=s" => \$file );

my $first_name;
my $second_name;
my $first_read;

$/ = '@';



open (INFILE, $file) or last;
$file =~ s/.fq//g;
my $outfile1 = $file."_R1.fq";
my $outfile2 = $file."_R2.fq";
open (OUTFILE1, ">$outfile1");
open (OUTFILE2, ">$outfile2");
while (my $read = <INFILE>)
	{
	chomp($read);
	unless ($read =~ m/.+/) {next;}
	my @line = split /\n/, $read;
	if ($line[0] =~ m/1$/) {($first_name = $line[0]) =~ s/1$//; $first_read = $read;}
	if ($line[0] =~ m/2$/) 
		{
		($second_name = $line[0]) =~ s/2$//; 
		my $second_read = $read;
		if ($second_name eq $first_name)
			{
			print OUTFILE1 "@".$first_read;
			print OUTFILE2 "@".$second_read;
			}
		$first_name= "";	
		}
	}	
