#!/usr/bin/perl

use File::Basename;
use Compress::Zlib;

$file=fileparse($ARGV[0]);
$outdir=$ARGV[1];
$normFacType=$ARGV[2];
$normFac=$ARGV[3];
if($file=~/gz/)
{
	$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
	while($gz->gzreadline($_) > 0)
	{
	chomp;
	split(/\t/);
	if($normFac)
	{
		$nf=$normFac/1000000;
		$sense{length($_[0])}+=$_[1]/$_[6]/$_[7]/$nf if ($_[3] eq 'sense');
		$antisense{length($_[0])}+=$_[1]/$_[6]/$_[7]/$nf if ($_[3] eq 'antisense');
	}
	else
	{
		$sense{length($_[0])}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
		$antisense{length($_[0])}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');
	}
	}
}
else
{
    open IN, $ARGV[0];
	while(<IN>)
	{
		chomp;
		split(/\t/);
		if($normFac)
		{
			$nf=$normFac/1000000;
			$sense{length($_[0])}+=$_[1]/$_[6]/$_[7]/$nf if ($_[3] eq 'sense');
			$antisense{length($_[0])}+=$_[1]/$_[6]/$_[7]/$nf if ($_[3] eq 'antisense');
		}
		else
		{
			$sense{length($_[0])}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
			$antisense{length($_[0])}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');
		}
	}
}

$file=~/(.*).gz/;
$filename=$1;
if($normFac)
{
	open OUT, ">$outdir/$filename.$normFacType.lendis2";
}
else
{
	open OUT, ">$outdir/$filename.lendis2";
}
foreach (18..29)
 {
 $sense{$_}=0 if (!exists $sense{$_});
 $antisense{$_}=0 if (!exists $antisense{$_});
 print OUT "$_\t$sense{$_}\t$antisense{$_}\n";
 }
 close(OUT);