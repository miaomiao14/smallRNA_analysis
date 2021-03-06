#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;

	my $file=$ARGV[0];
	my $filename=fileparse($file);
	my $OUTDIR=$ARGV[1];
	my %filehash=();
	if($file=~/gz/)
	{
		my $gz="";
		$gz = gzopen($file, "rb") or die "Cannot open $file: $gzerrno\n" ;
		while($gz->gzreadline(my $record) > 0)
		{ 
			chomp $record;
			my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
			$filehash{$chr}{$record}=1;
		}
		$gz->gzclose();
	}
	else
	{
		open IN, $file or die "Cannot open $file: $!\n";
		while(my $record=<IN>) 
		{
			chomp $record;
			my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
			$filehash{$chr}{$record}=1;
		}
		close(IN);
	}
	open CHR, ">$OUTDIR/$filename.chrlist.txt";
	foreach my $chr (keys %filehash)
	{
		print CHR $chr, "\n";
		open OUT, ">$OUTDIR/$filename.$chr" or die "Cannot open $OUTDIR/$filename.$chr to write: $!\n";
		foreach my $record (keys %{$filehash{$chr}})
		{
			print OUT $record,"\n";
		}
		close(OUT);
	}
	close(CHR);