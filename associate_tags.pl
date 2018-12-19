#!/usr/bin/perl

##################
#
#052314 - Updated to account for multimapping oligos
#
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
	
		${$tags{$revcomp}}[0]++;
		${$tags{$revcomp}}[1] = -9 ;
		${$tags{$revcomp}}[2] = "-" ;
		${$tags{$revcomp}}[3] = "-" ;
		${$tags{$revcomp}}[4] = $inline[2];
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
my $cur_ori;
my $cur_flag;

my @tmp_id;
my @tmp_passflg;
my @tmp_ori;
my $collision_ok;
my $p;
my $m;
my $i;

while (<ENHANCERS>)
	{
	chomp;
	@inline = split("\t");
	
	$cur_tag = $inline[0];
	$cur_loc = $inline[1];
	$cur_flag = $inline[4];
	
	if($cur_flag > 0)
		{
		if(exists($tags{$cur_tag}))
			{
			if($cur_flag == 1)
				{
				$collision_ok=1;
				@tmp_id = split(/,/,$cur_loc);
				$m = $tmp_id[0];
				$m =~ s/_RC//;
				foreach $p (@tmp_id)
					{
					$p =~ s/_RC//;
					$collision_ok=0 unless(exists($approved_multiHit{$tmp_id[0]}{$p}));
					}
				if($collision_ok == 0)
					{				
					${$tags{$cur_tag}}[1] = -4;
					${$tags{$cur_tag}}[2] = $cur_loc;
					${$tags{$cur_tag}}[3] = -9
					}
				if($collision_ok == 1)
					{
					@tmp_ori = ();
					for($i=0;$i<scalar(@tmp_id);$i++)
						{
						if($tmp_id[$i] =~ /_RC_/)
							{
							$tmp_id[$i] =~ s/_RC//;
							push(@tmp_ori,1);
							}
						else
							{
							push(@tmp_ori,0);							
							}
						}
					${$tags{$cur_tag}}[1] = -1;
					${$tags{$cur_tag}}[2] = join(",",@tmp_id);	
					${$tags{$cur_tag}}[3] = join(",",@tmp_ori)

					}	
				}
			elsif($cur_flag == 2)
				{
				${$tags{$cur_tag}}[1] = -5;
				${$tags{$cur_tag}}[2] = $cur_loc;
				${$tags{$cur_tag}}[3] = -9
				}
			elsif($cur_flag > 2)
				{
				${$tags{$cur_tag}}[1] = -6;
				${$tags{$cur_tag}}[2] = $cur_loc;
				${$tags{$cur_tag}}[3] = -9
				}
			}
		}
	else
		{
		if(exists($tags{$cur_tag}))
			{
			$cur_ori = 0;
			if($cur_loc =~ /_RC_/)
				{
				$cur_loc =~ s/_RC//;
				$cur_ori = 1;
				}
			${$tags{$cur_tag}}[1] = $cur_ori;
			${$tags{$cur_tag}}[2] = $cur_loc;
			${$tags{$cur_tag}}[3] = $cur_ori;		
			}
		}
	}	
close ENHANCERS;


my $key;

foreach $key (keys %tags)
	{
	print join("\t",$key,@{$tags{$key}})."\n";
	}
