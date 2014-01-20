#!/usr/bin/perl


use File::Basename;
use Compress::Zlib;
#open IN, $ARGV[0];



$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
while($gz->gzreadline($_) > 0) 
{ chomp; split(/\t/);
$_[2]=~/(\w.+):(\d+)-(\d+)\((.+)\)/;
$key="$1\t$2\t$3\t$4\t$_[0]\t$_[1]\t$_[6]\t$_[3]\n";
$hash{$key}=0;
}
$gz->gzclose();

$file=fileparse($ARGV[0]);
$file=~/(.*).gz/;
$filename=$1;
$OUTDIR=$ARGV[1];
open OUT,"> $OUTDIR/$filename.norm.bed";
foreach $record (keys %hash) { 
print OUT $record; 
}
close(OUT);
