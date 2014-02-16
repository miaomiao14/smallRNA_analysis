#!/usr/bin/perl

#Author: Jia Xu
#select a range of length of the seqs 
#lenselector <file> <len min> <len max>
use Compress::Zlib;
if(@ARGV<1) { print "USAGE:gzlenrangeselector.pl <file> <len1> <len2> \n"; exit 1;}
#open IN, $ARGV[0];
$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
while($gz->gzreadline($_) > 0)
{
 $line=$_;
 chomp $line;
 @l=split(/\t/,$line);
 print "$line\n" if (length($line[0])>=$ARGV[1] & length($line[0])<=$ARGV[2]);
} 
$gz->gzclose(); 