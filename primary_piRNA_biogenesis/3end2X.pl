#!/usr/bin/perl


##Find the nearest U from the 3' end of piRNA species###
##Record the length distribution of those walking distance###
##As a control, find the nearset A from the 3' end of piRNA species###

BEGIN { unshift @INC,"/home/xuj1/bin";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;



open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while (<IN>)
{
    if(/>(.+) type/)
    {
        $chr="chr$1";
    }
    else
    {
        chomp;
        $genome{$chr}=$_;
    }
}
close(IN);
open IN, "/home/wangw1/pipeline/common/dm3.chrom.sizes";
while(<IN>) {
    next if(/#/);
    chomp;split(/\t/);
             $chr=$_[0];   $chrsize{$chr}=$_[1];}
close(IN);
open IN, "$ARGV[0]";
open OUT, ">$ARGV[0].10ntaway.out";
%step=();
 while(<IN>) #complementarity level
{
    chomp;
    split(/\t/);
    next if (length($_[4])>29 || length($_[4])<23);
    #$total{$file}+=$_[5]/$_[6];
    next if (/data/);
   
    for($i=1;$i<=10;$i++)
    {
        if($_[3] eq '+')
        {
            $start=$_[2]+$i-1;
            $str=substr($genome{$_[0]},$start,1);
            #if($str eq $ARGV[1])
            #{
            $step{$str}{$i}+=$_[5]/$_[6];
            #last;
            #}
        }
        else
        {
            $start=$_[1]-$i+1; #for bed format file, the end is the actural end coordinate+1 see UCSC bed format
            $str=substr($genome{$_[0]},$start,1);
            $str=&revfa($str);
            #if($str eq $ARGV[1])
            #{
            $step{$str}{$i}+=$_[5]/$_[6];
            #last;
            #}
        }
    }
}


foreach $str (keys %step)
{
    foreach $distance (sort {$a <=> $b } keys %{$step{$str}} )
    {
        print OUT "$str\t$distance\t$step{$str}{$distance}\n";
    }
}
close(OUT);
close(IN);
