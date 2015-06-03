#!/usr/bin/perl


use File::Basename;
use Compress::Zlib;
#open IN, $ARGV[0];



$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
while($gz->gzreadline($line) > 0) 
{ chomp $line; @l=split(/\t/,$line);
$l[2]=~/(\w.+):(\d+)-(\d+)\((.+)\)/;
$key="$1\t$2\t$3\t$4\t$l[0]\t$l[1]\t$l[6]\t$l[3]\n";
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
