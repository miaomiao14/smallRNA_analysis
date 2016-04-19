#!/usr/bin/perl

BEGIN { unshift @INC,"/diag/home/netashawang/git/smallRNA_analysis/utils";}
require "restrict_digts.pm";
use File::Basename;

$bed2file=$ARGV[0];
$outfile=$ARGV[1];

my %mirstat;
my %mirtotal;
my %mirprecursor;

open BED, "$bed2file" or die "can not find file $bed2file: $!";
while(my $line=<BED>)
{
	chomp $line;
	my ($mir, $foffset, $toffset,$read,$ntm,$strand,$seq) = split ( /\t/, $line );
	

	$mirstat{$mir}{$foffset}+=$read/$ntm;
	$mirtotal{$mir}+=$read/$ntm;
	$mir =~ s/\-[53]p//;
	$mirprecursor{$mir}=1;
}
close(BED);

my %mirpercent;

foreach my $mir (keys %mirstat)
{
	foreach my $offset (keys %{$mirstat{$mir}})
	{
	$mirpercent{$mir}{$offset}=$mirstat{$mir}{$offset}/$mirtotal{$mir};
	$mirpercent{$mir}{$offset}=&restrict_num_decimal_digits($mirpercent{$mir}{$offset},5)
	}
}
open OUT, ">$outfile" or die "can not find file $outfile: $!";

foreach my $mirpre (keys %mirprecursor) ##redundancy
{

	#$mirmature =~ s/\-[53]p//;
	#$mirpre=$mirmature;
	#if($mirmature=~/5p/)
	my $mir5pname=$mirpre."-5p";
	my $mir3pname=$mirpre."-3p";
	$mir5pstat=$mirstat{$mir5pname}{0};
	$mir3pstat=$mirstat{$mir3pname}{0};
	my $gp=1;
	
	for (my $i=-1; $i<2; $i++)
	{
		if( !exists $mirpercent{$mir5pname}{$i})
		{
			$mirpercent{$mir5pname}{$i}=0;
		}
		if( !exists $mirpercent{$mir3pname}{$i})
		{
			$mirpercent{$mir3pname}{$i}=0;
		}
		
	}
    if($mir5pstat>$mir3pstat)
	{
		if($mir3pstat!=0)
		{
			$gp=$mir5pstat/$mir3pstat ;
			$gp=&restrict_num_decimal_digits($gp,3);
		}
		else
		{
			$gp="N/A";
		}

		
		print OUT $mir5pname,"\t",$gp,"\t","N-1:",$mirpercent{$mir5pname}{-1},"\t","N:",$mirpercent{$mir5pname}{0},"\t","N+1:",$mirpercent{$mir5pname}{1},"\n";
		
	}
	else
	{

		if($mir5pstat!=0)
		{
			$gp=$mir3pstat/$mir5pstat ;
			$gp=&restrict_num_decimal_digits($gp,3);
		}
		else
		{
			$gp="N/A";
		}
		
		print OUT $mir3pname,"\t",$gp,"\t","N-1:",$mirpercent{$mir3pname}{-1},"\t","N:",$mirpercent{$mir3pname}{0},"\t","N+1:",$mirpercent{$mir3pname}{1},"\n";
	}

}
close(OUT);