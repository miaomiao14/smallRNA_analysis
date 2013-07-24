#!/usr/bin/env perl
##This is a perl script trying to combine all the information from transposon mapper to circos data track
##input: DIR where the files are
##type mean or max or min
##binsize 
##prefix of file name
use File::Basename;
#my ($filename,$dir)=fileparse($ARGV[1]);
$dir=$ARGV[0];
$type=$ARGV[1];
@files=`ls $dir/*.$type.*`;
$binsize=$ARGV[2];
$pre=$ARGV[3];
@pre=split(/\./,$pre);
$pre=$pre[0];
open CHR, "/home/wangw1/pipeline/common/dm3.chrom.sizes";

while($read=<CHR>)
{
chomp $read;
@a=split(/\t/,$read);
$chrsize{$a[0]}=$a[1];
}
close(CHR);

open OUT, ">$dir/$pre.$type.circos.$binsize.txt";

foreach $f (@files)
{
open SCORE, "$f";
$sco=<SCORE>;
chomp $sco;
@scores=split(/\t/,$sco);
my ($filename,$dir)=fileparse($f);
print "$filename\n";
$filename=~/A.*.(chr.+)\.$type\.bin\.txt/;
#$r[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
$chrom=$1;
$sizechr=$chrsize{$chrom};
$nbins=$sizechr/$binsize;
	for($i=0;$i<($nbins-1);$i++)
	{
  	$start=$binsize*$i;
  	$end=$binsize*($i+1)-1;
 		if($scores[$i] && $scores[$i] ne "n/a")
  		{
  		$binscore=$scores[$i];
  		}
  		else
  		{
  		$binscore=0;
  		}
  	print OUT "$chrom\t$start\t$end\t$binscore\n";
	}
	if($i!=0)
	{
	$start=$binsize*($i);
	$end=$sizechr;
		if($scores[$i] && $scores[$i] ne "n/a")
  		{
  		$binscore=$scores[$i];
  		}
  		else
  		{
  		$binscore=0;
		}
  		print OUT "$chrom\t$start\t$end\t$binscore\n";
	}
	else
	{
	$start=0;
	$end=$sizechr;
		if($scores[$i] && $scores[$i] ne "n/a")
  		{
  		$binscore=$scores[$i];
  		}
  		else
  		{
  		$binscore=0;
  		}
  	print OUT "$chrom\t$start\t$end\t$binscore\n";   
	}
}
close(SCORE);
close(OUT);