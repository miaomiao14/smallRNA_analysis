#!/usr/bin/perl

use File::Basename;
open IN, $ARGV[0];
$filename=fileparse($ARGV[0]);
$OUTDIR=$ARGV[1];
while(<IN>) {
  chomp;
  split(/\t/);
  $_[5] =~ /(.+)\./;
  $transposon{$1}.="$_\n";
}
close(IN);

foreach $F(keys %transposon)
 {
  open OUT,"> $OUTDIR/$filename.$F";
  print OUT "$transposon{$F}";
  close(OUT);
}