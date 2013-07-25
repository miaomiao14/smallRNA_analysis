#!/usr/bin/perl

#Author: Jia Xu
#select a range of length of the seqs 
#lenselector <file> <len min> <len max>
use Compress::Zlib;
if(@ARGV<1) { print "USAGE:lenselector <file> <len1> <len2> \n"; exit 1;}
#open IN, $ARGV[0];
$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
while($gz->gzreadline($_) > 0)
{
 chomp;
 split(/\t/);
 print "$_\n" if (length($_[0])>=$ARGV[1] & length($_[0])<=$ARGV[2]);
} 