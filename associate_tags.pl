#!/usr/bin/perl

##################
#
#052314 - Updated to account for multimapping oligos
#
####
#	FLAGS
##################

use strict;
use warnings;
 
 
my $tags = $ARGV[0]; #matched tag file
my $enhancers = $ARGV[1];
my $multi_hit = $ARGV[2];

open (ENHANCERS, "$enhancers") or die("ERROR: can not read file ($enhancers): $!\n");
open (MULTIHIT, "$multi_hit") or die("ERROR: can not read file ($multi_hit): $!\n");

#open (OUT, ">$out") or die("ERROR: can not create $out: $!\n");

my @inline;
my %tags; #[ct,orientation,location,tagseq(rc)]
my $revcomp;

if($ARGV[0] eq "stdin")
{
while (<STDIN>)
	{
	
	chomp;
	@inline = split("\t");
	
	if($inline[6] eq "PASS")
		{
		$revcomp = reverse($inline[2]);
		$revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
	
		${$tags{$revcomp}}[0]++;  #Coverage
		${$tags{$revcomp}}[1] = -9 ; #Depreciated flag (for ori)
		${$tags{$revcomp}}[2] = "-" ; #Oligo ID
		${$tags{$revcomp}}[3] = "-" ;	#mapping flag
		${$tags{$revcomp}}[4] = $inline[2]; #RC Tag
		${$tags{$revcomp}}[5] = "NA"; #mapping score
		${$tags{$revcomp}}[6] = "NA"; #CIGAR (future)
		${$tags{$revcomp}}[7] = "NA"; #MD tag (future)


		}
	}
}
else
{
open (TAGS, "$tags") or die("ERROR: can not read file ($tags): $!\n");
while (<TAGS>)
	{
	
	chomp;
	@inline = split("\t");
	
	if($inline[6] eq "PASS")
		{
		$revcomp = reverse($inline[2]);
		$revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
	
		${$tags{$revcomp}}[0]++;
		${$tags{$revcomp}}[1] = -9 ;
		${$tags{$revcomp}}[2] = "-" ;
		${$tags{$revcomp}}[3] = "-" ;
		${$tags{$revcomp}}[4] = $inline[2];
		${$tags{$revcomp}}[5] = "NA"; #mapping score
		${$tags{$revcomp}}[6] = "NA"; #CIGAR (future)
		${$tags{$revcomp}}[7] = "NA"; #MD tag (future)
		}
	}
close TAGS;
}

my %approved_multiHit;
my $id;

while (<MULTIHIT>)
	{
	chomp;
	@inline = split("\t");
	foreach $id (@inline)
		{
		$approved_multiHit{$inline[0]}{$id}=1;
		}
	}


my $cur_tag;
my $cur_loc;
my $cur_flag;
my $cur_m_flag;
my $cur_m_aln;
my $cur_m_cigar;
my $cur_m_md;
	
my @tmp_id;
my @tmp_passflg;
my @tmp_aln;
my @tmp_cigar;
my @tmp_md;
my @tmp_ori;
my $collision_ok;
my $p;
my $i;

while (<ENHANCERS>)
	{
	chomp;
	@inline = split("\t");
	
	$cur_tag = $inline[0];
	$cur_loc = $inline[1];
	$cur_flag = $inline[4];
	$cur_m_flag = $inline[5];
	$cur_m_aln = $inline[6];
	$cur_m_cigar = $inline[7];
	$cur_m_md = $inline[8];
	
	if($cur_flag > 0)
		{
		if(exists($tags{$cur_tag}))
			{
			if($cur_flag == 1)
				{
				$collision_ok=1;
				@tmp_id = split(/,/,$cur_loc);
				@tmp_passflg = split(/,/,$cur_m_flag);
				@tmp_aln = split(/,/,$cur_m_aln);
				@tmp_cigar = split(/,/,$cur_m_cigar);
				@tmp_md = split(/,/,$cur_m_md);
				
				foreach $p (@tmp_id)
					{
					$collision_ok=0 unless(exists($approved_multiHit{$tmp_id[0]}{$p}));
					}
				if($collision_ok == 0)
					{				
					${$tags{$cur_tag}}[1] = -4;
					${$tags{$cur_tag}}[2] = $cur_loc;
					${$tags{$cur_tag}}[3] = $cur_flag;
					${$tags{$cur_tag}}[5] = $cur_m_aln;
					${$tags{$cur_tag}}[6] = $cur_m_cigar; #CIGAR
					${$tags{$cur_tag}}[7] = $cur_m_md; #MD tag
					}
				if($collision_ok == 1)
					{
					for($i=0;$i<scalar(@tmp_id);$i++)
						{
						${$tags{$cur_tag}}[1] = -1;
						${$tags{$cur_tag}}[2] = $tmp_id[$i];	
						${$tags{$cur_tag}}[3] = $tmp_passflg[$i];
						${$tags{$cur_tag}}[5] = $tmp_aln[$i];
						${$tags{$cur_tag}}[6] = $tmp_cigar[$i]; #CIGAR
						${$tags{$cur_tag}}[7] = $tmp_md[$i]; #MD tag
						}
					}	
				}
			elsif($cur_flag == 2)
				{
				${$tags{$cur_tag}}[1] = -5;
				${$tags{$cur_tag}}[2] = $cur_loc;
				${$tags{$cur_tag}}[3] = $cur_flag;
				${$tags{$cur_tag}}[5] = $cur_m_aln;
				${$tags{$cur_tag}}[6] = $cur_m_cigar; #CIGAR
				${$tags{$cur_tag}}[7] = $cur_m_md; #MD tag
				}
			elsif($cur_flag > 2)
				{
				${$tags{$cur_tag}}[1] = -6;
				${$tags{$cur_tag}}[2] = $cur_loc;
				${$tags{$cur_tag}}[3] = $cur_flag;
				${$tags{$cur_tag}}[5] = $cur_m_aln;
				${$tags{$cur_tag}}[6] = $cur_m_cigar; #CIGAR
				${$tags{$cur_tag}}[7] = $cur_m_md; #MD tag
				}
			}
		}
	else
		{
		if(exists($tags{$cur_tag}))
			{
			${$tags{$cur_tag}}[1] = 0;
			${$tags{$cur_tag}}[2] = $cur_loc;
			${$tags{$cur_tag}}[3] = $cur_flag;	
			${$tags{$cur_tag}}[5] = $cur_m_aln;
			${$tags{$cur_tag}}[6] = $cur_m_cigar;
			${$tags{$cur_tag}}[7] = $cur_m_md;	
			}
		}
	}	
close ENHANCERS;


my $key;

foreach $key (keys %tags)
	{
	print join("\t",$key,@{$tags{$key}})."\n";
	}
