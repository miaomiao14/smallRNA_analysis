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
	
	$rooReads=0;
	%roo=();
	while($gz->gzreadline($_) > 0)
	{
	chomp;
	split(/\t/);
	#needs to calculate the reads mapped to roo, update the nf;
	if(/FBgn0000155_roo/){$roo{$_[0]}=$_[1];}
	}
	foreach $species (keys %roo)
	{
		$rooReads+=$roo{$species};
	} 
	
	
	while($gz->gzreadline($_) > 0)
	{
	chomp;
	split(/\t/);
	next if(/FBgn0000155_roo/);

	if($normFac)
	{
		$nf=($normFac-$rooReads)/1000000;
		$sense{length($_[0])}+=$_[1]/$_[6]/$_[7]/$nf if ($_[3] eq 'sense');
		$antisense{length($_[0])}+=$_[1]/$_[6]/$_[7]/$nf if ($_[3] eq 'antisense');
	}
	else
	{
		$sense{length($_[0])}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
		$antisense{length($_[0])}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');
	}
	}
	$gz->gzclose();
}
else
{
    open IN, $ARGV[0];
    $rooReads=0;
	%roo=();
    while(<IN>)
	{
		chomp;
		split(/\t/);
		if(/FBgn0000155_roo/){$roo{$_[0]}=$_[1];}
	}
    foreach $species (keys %roo)
	{
		$rooReads+=$roo{$species};
	}
    
	while(<IN>)
	{
		chomp;
		split(/\t/);
		next if(/FBgn0000155_roo/);
		if($normFac)
		{
			$nf=($normFac-$rooReads)/1000000;
			$sense{length($_[0])}+=$_[1]/$_[6]/$_[7]/$nf if ($_[3] eq 'sense');
			$antisense{length($_[0])}+=$_[1]/$_[6]/$_[7]/$nf if ($_[3] eq 'antisense');
		}
		else
		{
			$sense{length($_[0])}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
			$antisense{length($_[0])}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');
		}
	}
	close(IN);
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