#!/usr/bin/perl

use File::Basename;
use Compress::Zlib;
#open IN, $ARGV[0];
$file=fileparse($ARGV[0]);
$file=~/(.*).gz/;
$filename=$1;
$OUTDIR=$ARGV[1];
$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
while($gz->gzreadline($_) > 0)
#while(<IN>) 
{
  chomp;
  split(/\t/);
  $_[5] =~ /(.+)\./;
  $transposon{$1}.="$_\n";
}
#close(IN);
$gz->gzclose();
foreach $F(keys %transposon)
 {
  open OUT,"> $OUTDIR/$filename.$F";
  print OUT "$transposon{$F}";
  close(OUT);
}