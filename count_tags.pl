#!/usr/bin/perl

use strict;
use warnings;
use Statistics::Basic qw(:all);

 
my $list = $ARGV[0]; #matched tag file
my $counts = $ARGV[1];
my $match_file = $ARGV[2];
my $cutoff = $ARGV[3];
my $out = $ARGV[4];

open (LIST, "$list") or die("ERROR: can not read file ($list): $!\n");
open (MATCH, "$match_file") or die("ERROR: can not read file ($match_file): $!\n");
open (COUNTS, "$counts") or die("ERROR: can not read file ($counts): $!\n");

#open (OUT, ">$out") or die("ERROR: can not create $out: $!\n");

my @inline;
my %pos_strand_ct; 
my %neg_strand_ct;
my %pos_strand_indiv;
my %neg_strand_indiv;
my %pos_strand_indiv_tag;
my %neg_strand_indiv_tag;
my %pos_strand_events;
my %neg_strand_events;
my %order;
my $p = 1;

while (<LIST>)
	{
	chomp;
	@inline = split("\t");
	$pos_strand_ct{$inline[0]}=0;
	$neg_strand_ct{$inline[0]}=0;
	@{$pos_strand_indiv{$inline[0]}}= ();
	@{$neg_strand_indiv{$inline[0]}}= ();
	@{$pos_strand_indiv_tag{$inline[0]}}= ();
	@{$neg_strand_indiv_tag{$inline[0]}}= ();
	$pos_strand_events{$inline[0]}=0;
	$neg_strand_events{$inline[0]}=0;
	$order{$inline[0]}=$p;
	$p++;
	}
close LIST;

my %matches;

while (<MATCH>)
	{
	chomp;
	@inline = split("\t");
	$matches{$inline[0]}=$inline[1];	
	}

my $tag;
my $ct;
my $strand;
my $id;
my $events = 0;
my $tag_events = 0;
my $tag_events_greater1 = 0;
my $total;
my $total_all = 0;
my $total_all_greater1 = 0;
my @multi_id;
my @multi_strand;
my %multi_id_hash;
my $id_m;
my $i_ct;
my $i_first;

while (<COUNTS>)
	{
	chomp;
	@inline = split("\t");
	#$tag = $inline[0];
	$tag = ConvertDNAtoNUM($inline[0]);
	$ct = $inline[1];
	$strand = $inline[2];
	$id = $inline[3];
	$tag_events++;
	$tag_events_greater1++ if($ct > 1);
	$total_all += $ct;
	$total_all_greater1 += $ct if($ct > 1);
	
	if($ct >= $cutoff)
		{
		if($strand == 0)
			{
			$pos_strand_ct{$id} += $ct;
			push(@{$pos_strand_indiv{$id}},$ct);
			push(@{$pos_strand_indiv_tag{$id}},$tag);
			$pos_strand_events{$id}++;
			$events++;
			$total+=$ct;
			if(exists $matches{$id})
				{
				$pos_strand_ct{$matches{$id}} += $ct;
				push(@{$pos_strand_indiv{$matches{$id}}},$ct);
				push(@{$pos_strand_indiv_tag{$matches{$id}}},$tag);
				$pos_strand_events{$matches{$id}}++;
				}
			}
		elsif($strand == 1)
			{
			$neg_strand_ct{$id} += $ct;	
			push(@{$neg_strand_indiv{$id}},$ct);
			push(@{$neg_strand_indiv_tag{$id}},$tag);
			$neg_strand_events{$id}++;
			$events++;
			$total+=$ct;
			if(exists $matches{$id})
				{
				$neg_strand_ct{$matches{$id}} += $ct;
				push(@{$neg_strand_indiv{$matches{$id}}},$ct);
				push(@{$neg_strand_indiv_tag{$matches{$id}}},$tag);
				$neg_strand_events{$matches{$id}}++;
				}
			}
		elsif($strand == -1)
			{
			@multi_id = split(/,/,$inline[3]);
			@multi_strand = split(/,/,$inline[4]);
			foreach $id_m (@multi_id){$multi_id_hash{$id_m}=1;}
			$i_first=0;
			$i_ct=0;
			foreach $id_m (@multi_id)
				{
				$i_first++;
				if($multi_strand[$i_ct] == 0)
					{
					$pos_strand_ct{$id_m} += $ct;
					push(@{$pos_strand_indiv{$id_m}},$ct);
					push(@{$pos_strand_indiv_tag{$id_m}},$tag);
					$pos_strand_events{$id_m}++;
					$events++ if($i_first == 1);
					$total+=$ct if($i_first == 1);
					if(exists $matches{$id_m})
						{
						unless(exists $multi_id_hash{$matches{$id_m}})
							{
							$pos_strand_ct{$matches{$id_m}} += $ct;
							push(@{$pos_strand_indiv{$matches{$id_m}}},$ct);
							push(@{$pos_strand_indiv_tag{$matches{$id_m}}},$tag);
							$pos_strand_events{$matches{$id_m}}++;
							}
						}
					}
				elsif($multi_strand[$i_ct] == 1)
					{
					$neg_strand_ct{$id_m} += $ct;	
					push(@{$neg_strand_indiv{$id_m}},$ct);
					push(@{$neg_strand_indiv_tag{$id_m}},$tag);
					$neg_strand_events{$id_m}++;
					$events++ if($i_first == 1);
					$total+=$ct if($i_first == 1);
					if(exists $matches{$id_m})
						{
						unless(exists $multi_id_hash{$matches{$id_m}})
							{
							$neg_strand_ct{$matches{$id_m}} += $ct;
							push(@{$neg_strand_indiv{$matches{$id_m}}},$ct);
							push(@{$neg_strand_indiv_tag{$matches{$id_m}}},$tag);
							$neg_strand_events{$matches{$id_m}}++;
							}
						}
					}
				$i_ct++;
				}			
			}
		else
			{
		
			}
		}

	}	

my $pos_fpkm;
my $neg_fpkm;
my $total_fpkm;
#my $pos_event_fpkm_mean;
#my $neg_event_fpkm_mean;
#my $total_event_fpkm_mean;
#my $pos_event_fpkm_med;
#my $neg_event_fpkm_med;
#my $total_event_fpkm_med;
my $tmp_ct;
my $tmp_fpkm;
my %pos_strand_indiv_fpkm;
my %neg_strand_indiv_fpkm;


open (OUT_raw, ">$out.comb.raw") or die("ERROR: can not create $out: $!\n");
open (OUT_rpm, ">$out.comb.rpm") or die("ERROR: can not create $out: $!\n");
open (OUT_rpm_all, ">$out.comb.all.rpm") or die("ERROR: can not create $out: $!\n");
open (OUT_rpm_one, ">$out.comb.one.rpm") or die("ERROR: can not create $out: $!\n");
open (OUT_tpmm, ">$out.comb.tpmm") or die("ERROR: can not create $out: $!\n");
open (OUT_tpmm_all, ">$out.comb.all.tpmm") or die("ERROR: can not create $out: $!\n");
open (OUT_tpmm_one, ">$out.comb.one.tpmm") or die("ERROR: can not create $out: $!\n");

open (OUT_tag_id, ">$out.indiv.tagID") or die("ERROR: can not create $out: $!\n");
open (OUT_raw_tag, ">$out.indiv.raw") or die("ERROR: can not create $out: $!\n");
open (OUT_rpm_tag, ">$out.indiv.rpm") or die("ERROR: can not create $out: $!\n");
open (OUT_rpm_all_tag, ">$out.indiv.all.rpm") or die("ERROR: can not create $out: $!\n");
open (OUT_rpm_one_tag, ">$out.indiv.one.rpm") or die("ERROR: can not create $out: $!\n");
open (OUT_tpmm_tag, ">$out.indiv.tpmm") or die("ERROR: can not create $out: $!\n");
open (OUT_tpmm_all_tag, ">$out.indiv.all.tpmm") or die("ERROR: can not create $out: $!\n");
open (OUT_tpmm_one_tag, ">$out.indiv.one.tpmm") or die("ERROR: can not create $out: $!\n");
#######
##Print Raw Counts
#[+ strand ct,- strand ct,combind ct]
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
print OUT_raw join ("\t",$pos_strand_ct{$id},$neg_strand_ct{$id},$pos_strand_ct{$id}+$neg_strand_ct{$id})."\n";
}


##Print RPM corrected for depth at associated tags
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
$pos_fpkm = sprintf("%.5f", ($pos_strand_ct{$id}/$total)*1000000);
$neg_fpkm = sprintf("%.5f", ($neg_strand_ct{$id}/$total)*1000000);
$total_fpkm = sprintf("%.5f", (($neg_strand_ct{$id}+$pos_strand_ct{$id})/$total)*1000000);
#[,+ strand FPKM, - strand FPKM, combine FPKM]
print OUT_rpm join ("\t",$pos_fpkm,$neg_fpkm,$total_fpkm)."\n";
}

##Print RPM corrected for depth at all tags
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
$pos_fpkm = sprintf("%.5f", ($pos_strand_ct{$id}/$total_all)*1000000);
$neg_fpkm = sprintf("%.5f", ($neg_strand_ct{$id}/$total_all)*1000000);
$total_fpkm = sprintf("%.5f", (($neg_strand_ct{$id}+$pos_strand_ct{$id})/$total_all)*1000000);
#[,+ strand FPKM, - strand FPKM, combine FPKM]
print OUT_rpm_all join ("\t",$pos_fpkm,$neg_fpkm,$total_fpkm)."\n";
}

##Print RPM corrected for depth at tags with greater than one count
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
$pos_fpkm = sprintf("%.5f", ($pos_strand_ct{$id}/$total_all_greater1)*1000000);
$neg_fpkm = sprintf("%.5f", ($neg_strand_ct{$id}/$total_all_greater1)*1000000);
$total_fpkm = sprintf("%.5f", (($neg_strand_ct{$id}+$pos_strand_ct{$id})/$total_all_greater1)*1000000);
#[,+ strand FPKM, - strand FPKM, combine FPKM]
print OUT_rpm_one join ("\t",$pos_fpkm,$neg_fpkm,$total_fpkm)."\n";
}

##Print TPMM (reads & tag) corrected for depth for depth at associated tags
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
$pos_fpkm = sprintf("%.5f", ($pos_strand_ct{$id}*(1000000/$total)*(1000000/$events)));
$neg_fpkm = sprintf("%.5f", ($neg_strand_ct{$id}*(1000000/$total)*(1000000/$events)));
$total_fpkm = sprintf("%.5f", (($neg_strand_ct{$id}+$pos_strand_ct{$id})*(1000000/$total)*(1000000/$events)));
#[,+ strand FPKM, - strand FPKM, combine FPKM]
print OUT_tpmm join ("\t",$pos_fpkm,$neg_fpkm,$total_fpkm)."\n";
}

##Print TPMM (reads & tag) corrected for depth at all tags
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
$pos_fpkm = sprintf("%.5f", ($pos_strand_ct{$id}*(1000000/$total_all)*(1000000/$tag_events)));
$neg_fpkm = sprintf("%.5f", ($neg_strand_ct{$id}*(1000000/$total_all)*(1000000/$tag_events)));
$total_fpkm = sprintf("%.5f", (($neg_strand_ct{$id}+$pos_strand_ct{$id})*(1000000/$total_all)*(1000000/$tag_events)));
#[,+ strand FPKM, - strand FPKM, combine FPKM]
print OUT_tpmm_all join ("\t",$pos_fpkm,$neg_fpkm,$total_fpkm)."\n";
}

##Print TPMM (reads & tag) corrected for depth at all tags with greater than one count
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
$pos_fpkm = sprintf("%.5f", ($pos_strand_ct{$id}*(1000000/$total_all_greater1)*(1000000/$tag_events_greater1)));
$neg_fpkm = sprintf("%.5f", ($neg_strand_ct{$id}*(1000000/$total_all_greater1)*(1000000/$tag_events_greater1)));
$total_fpkm = sprintf("%.5f", (($neg_strand_ct{$id}+$pos_strand_ct{$id})*(1000000/$total_all_greater1)*(1000000/$tag_events_greater1)));
#[,+ strand FPKM, - strand FPKM, combine FPKM]
print OUT_tpmm_one join ("\t",$pos_fpkm,$neg_fpkm,$total_fpkm)."\n";
}


#######
##Print individual tag counts for pos and neg strands
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
print OUT_raw_tag join ("\t",join(",","@",@{$pos_strand_indiv{$id}}),join(",","@",@{$neg_strand_indiv{$id}}))."\n";
print OUT_tag_id join ("\t",join(",","@",@{$pos_strand_indiv_tag{$id}}),join(",","@",@{$neg_strand_indiv_tag{$id}}))."\n";

}

##Print RPM corrected for depth at associated tags
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
@{$pos_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$pos_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct/$total*1000000));
push(@{$pos_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
@{$neg_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$neg_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct/$total*1000000));
push(@{$neg_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
print OUT_rpm_tag join ("\t",join(",","@",@{$pos_strand_indiv_fpkm{$id}}),join(",","@",@{$neg_strand_indiv_fpkm{$id}}))."\n";
}

##Print RPM corrected for depth at all tags
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
@{$pos_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$pos_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct/$total_all*1000000));
push(@{$pos_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
@{$neg_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$neg_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct/$total_all*1000000));
push(@{$neg_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
print OUT_rpm_all_tag join ("\t",join(",","@",@{$pos_strand_indiv_fpkm{$id}}),join(",","@",@{$neg_strand_indiv_fpkm{$id}}))."\n";
}

##Print RPM corrected for depth at tags with greater than one count
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
@{$pos_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$pos_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct/$total_all_greater1*1000000));
push(@{$pos_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
@{$neg_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$neg_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct/$total_all_greater1*1000000));
push(@{$neg_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
print OUT_rpm_one_tag join ("\t",join(",","@",@{$pos_strand_indiv_fpkm{$id}}),join(",","@",@{$neg_strand_indiv_fpkm{$id}}))."\n";
}

##Print TPMM (reads & tag) corrected for depth for depth at associated tags
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
@{$pos_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$pos_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct*(1000000/$total)*(1000000/$events)));
push(@{$pos_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
@{$neg_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$neg_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct*(1000000/$total)*(1000000/$events)));
push(@{$neg_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
print OUT_tpmm_tag join ("\t",join(",","@",@{$pos_strand_indiv_fpkm{$id}}),join(",","@",@{$neg_strand_indiv_fpkm{$id}}))."\n";
}

##Print TPMM (reads & tag) corrected for depth at all tags
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
@{$pos_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$pos_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct*(1000000/$total_all)*(1000000/$tag_events)));
push(@{$pos_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
@{$neg_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$neg_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct*(1000000/$total_all)*(1000000/$tag_events)));
push(@{$neg_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
print OUT_tpmm_all_tag join ("\t",join(",","@",@{$pos_strand_indiv_fpkm{$id}}),join(",","@",@{$neg_strand_indiv_fpkm{$id}}))."\n";
}

##Print TPMM (reads & tag) corrected for depth at all tags with greater than one count
foreach $id (sort { $order{$a} <=> $order{$b} } keys(%order))
{
@{$pos_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$pos_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct*(1000000/$total_all_greater1)*(1000000/$tag_events_greater1)));
push(@{$pos_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
@{$neg_strand_indiv_fpkm{$id}} = ();
foreach $tmp_ct (@{$neg_strand_indiv{$id}})
{
$tmp_fpkm = sprintf("%.5f", ($tmp_ct*(1000000/$total_all_greater1)*(1000000/$tag_events_greater1)));
push(@{$neg_strand_indiv_fpkm{$id}},$tmp_fpkm);
}
print OUT_tpmm_one_tag join ("\t",join(",","@",@{$pos_strand_indiv_fpkm{$id}}),join(",","@",@{$neg_strand_indiv_fpkm{$id}}))."\n";
}

#$pos_event_fpkm_mean = 0;
#$neg_event_fpkm_mean = 0;
#$total_event_fpkm_mean = 0;
#$pos_event_fpkm_mean = sprintf("%.2f",$pos_event_fpkm_mean/$pos_strand_events{$id}) if($pos_strand_events{$id} > 0);
#$neg_event_fpkm_mean = sprintf("%.2f",$neg_event_fpkm_mean/$neg_strand_events{$id}) if($neg_strand_events{$id} > 0);
#$total_event_fpkm_mean = sprintf("%.2f",$total_event_fpkm_mean/($pos_strand_events{$id}+$neg_strand_events{$id})) if(($pos_strand_events{$id}+$neg_strand_events{$id}) > 0);

#$pos_event_fpkm_med = median(@{$pos_strand_indiv_fpkm{$id}});
#$neg_event_fpkm_med = median(@{$neg_strand_indiv_fpkm{$id}});
#$total_event_fpkm_med = median((@{$pos_strand_indiv_fpkm{$id}},@{$neg_strand_indiv_fpkm{$id}}));

#[+ strand ct,- strand ct,combind ct,+ strand FPKM, - strand FPKM, combine FPKM,pos events, neg events, total events, all events + strand indiv counts, - strand indiv counts]
#print join("\t",$id)."\t";
#print join ("\t",$pos_strand_ct{$id},$neg_strand_ct{$id},$pos_strand_ct{$id}+$neg_strand_ct{$id})."\t";
#print join ("\t",$pos_fpkm,$neg_fpkm,$total_fpkm)."\t";
#print join ("\t",$pos_strand_events{$id},$neg_strand_events{$id},$pos_strand_events{$id}+$neg_strand_events{$id})."\t";
#print join ("\t",$pos_event_fpkm_mean,$neg_event_fpkm_mean,$total_event_fpkm_mean,$events)."\t";
#print join ("\t",$pos_event_fpkm_med,$neg_event_fpkm_med,$total_event_fpkm_med)."\t";
#print join ("\t",join(",","@",@{$pos_strand_indiv{$id}}),join(",","@",@{$neg_strand_indiv{$id}}))."\t";
#print join ("\t",join(",","@",@{$pos_strand_indiv_fpkm{$id}}),join(",","@",@{$neg_strand_indiv_fpkm{$id}}))."\n";

#[+ strand ct,- strand ct,combind ct,+ strand FPKM, - strand FPKM, combine FPKM, + strand indiv counts, - strand indiv counts]
#print join("\t",$id)."\t";
#print join ("\t",$pos_strand_ct{$id},$neg_strand_ct{$id},$pos_strand_ct{$id}+$neg_strand_ct{$id})."\t";
#print join ("\t",$pos_fpkm,$neg_fpkm,$total_fpkm)."\t";
#print join ("\t",join(",","@",@{$pos_strand_indiv{$id}}),join(",","@",@{$neg_strand_indiv{$id}}))."\t";
#print join ("\t",join(",","@",@{$pos_strand_indiv_fpkm{$id}}),join(",","@",@{$neg_strand_indiv_fpkm{$id}}))."\n";



sub ConvertDNAtoNUM{
	no warnings 'portable';
   # Private variable for PrintHello function
	my @seq=split(//,$_[0]);
	my $bin="";
	foreach my $base (@seq)
		{	
		if($base eq "A"){$bin=$bin."00";}
		if($base eq "C"){$bin=$bin."01";}
		if($base eq "G"){$bin=$bin."10";}
		if($base eq "T"){$bin=$bin."11";}
		}
	my $numeric=oct("0b".$bin);
	return($numeric)
}
