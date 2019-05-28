#!/usr/bin/perl

use strict;
use warnings;
use Text::LevenshteinXS qw(distance);
use Getopt::Std;


my %options=();
getopts('HA', \%options);

#####
#
#-H = No Adapter search. Do a hard trim
#-A = No anchor search. 
#
#####

my $skip_adapter;

if(exists($options{H}))
	{
	print STDERR "Skipping adapter search\n";
	$skip_adapter = 1;
	}
else {$skip_adapter = 0;}

my $skip_anchor;

if(exists($options{R}))
	{
	print STDERR "Skipping anchor search\n";
	$skip_anchor = 1;
	}
else {$skip_anchor = 0;}

 
my $BARCODE_SIZE = $ARGV[0];
my $L1_SEQ_FULL = $ARGV[1];
my $out = $ARGV[2];

#open (FASTA, "$fasta") or die("ERROR: can not read file ($fasta): $!\n");
open (MATCH, ">$out".".match") or die("ERROR: can not create $out .matched: $!\n");
open (REJECT, ">$out".".reject") or die("ERROR: can not create $out .rejected: $!\n");
			   
###
#
# v3 mpra design TCTAGAGGTTCGTCG
#
#
#
#			   
my $L1_SEQ = $L1_SEQ_FULL;
#my $BARCODE_SIZE = 20;
my $L1_SEQ_ANCHOR = substr($L1_SEQ,0,2);


my $r1;
my $id;
my $L1_read_seq;
my $dist;

my $barcode_seq;

my $match1 = -9;

#my $match_dist = 4; #max levenshtein distance for match
my $match_dist = 4; #max levenshtein distance for match


my $i;
my $read_Length;


while (<STDIN>)
{
	chomp;
	$id=$_;
	$id =~ s/^>//;
	$id =~ s/\/1$//;

	$r1=<STDIN>;
	chomp $r1;
        $r1=substr($r1,0,$BARCODE_SIZE+length($L1_SEQ_FULL));

##Right side adapter match
	$match1 = -9;
	$dist = -9;
	$L1_SEQ = substr($L1_SEQ_FULL,0,(length($r1)-$BARCODE_SIZE));
	$match_dist = length($L1_SEQ)*.30;
	
	if(substr($r1,$BARCODE_SIZE,length($L1_SEQ_ANCHOR)) eq $L1_SEQ_ANCHOR) #Anchor match
		{
		$match1 = 0;		
		if(abs(distance($L1_SEQ, substr($r1,$BARCODE_SIZE,length($L1_SEQ)))) <= $match_dist)  ##Adapter match
				{
				$dist = distance($L1_SEQ, substr($r1,$BARCODE_SIZE,length($L1_SEQ)));
				$match1 = 0;
				#Matched internal adapter
				}
			else
				{
				$dist = distance($L1_SEQ, substr($r1,$BARCODE_SIZE,length($L1_SEQ)));
				$match1 = -999;
				$match1 = 0 if($skip_anchor == 1);
				}
		}
			
	if($match1 != -9)
		{
		$barcode_seq = substr($r1,0,$BARCODE_SIZE);
		$L1_read_seq = substr($r1,$BARCODE_SIZE,length($L1_SEQ));

		if($match1 == 0 && $barcode_seq !~ /N/)
			{
			print MATCH join("\t",$id,$match1,$barcode_seq,$L1_read_seq,$dist,$r1,"PASS")."\n";
			}
		else
			{
			print REJECT join("\t",$id,$match1,$barcode_seq,$L1_read_seq,$dist,$r1,"REJECT")."\n";		
			}
		}
	elsif ($skip_anchor == 1)
		{
		$barcode_seq = substr($r1,0,$BARCODE_SIZE);
		$L1_read_seq = substr($r1,$BARCODE_SIZE,length($L1_SEQ));
		if($barcode_seq !~ /N/)
			{
			print MATCH join("\t",$id,$match1,$barcode_seq,$L1_read_seq,$dist,$r1,"PASS")."\n";
			}
		else
			{
			print REJECT join("\t",$id,$match1,$barcode_seq,$L1_read_seq,$dist,$r1,"REJECT")."\n";		
			}	
		}
	else
		{
		$barcode_seq = substr($r1,0,$BARCODE_SIZE);
		$L1_read_seq = substr($r1,$BARCODE_SIZE,length($L1_SEQ));
		print REJECT join("\t",$id,$match1,$barcode_seq,$L1_read_seq,$dist,$r1,"REJECT")."\n";	
		}
$r1=<STDIN>;
$r1=<STDIN>;		
}

