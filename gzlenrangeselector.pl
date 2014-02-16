#!/usr/bin/perl

#Author: Jia Xu
#select a range of length of the seqs 
#lenselector <file> <len min> <len max>
use Compress::Zlib;
if(@ARGV<1) { print "USAGE:gzlenrangeselector.pl <file> <len1> <len2> \n"; exit 1;}
#open IN, $ARGV[0];
$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
$infile=$ARGV[0];
$infile=~s/.gz//;
$outfile=$infile.".".$ARGV[1];
$outfile=$outfile."-".$ARGV[2]."gz";
$gzw=gzopen($outfile,"wb") or die "Cannot open $outfile: $gzerrno\n";
while($gz->gzreadline($line) > 0)
{
 #$line=$_;
 chomp $line;
 @l=split(/\t/,$line);
 $gzw->gzwrite( "$line\n") if (length($l[0])>=$ARGV[1] & length($l[0])<=$ARGV[2]);
} 
$gz->gzclose();
$gzw->gzclose();  