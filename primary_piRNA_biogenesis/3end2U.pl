#!/usr/bin/perl


##Find the nearest U from the 3' end of piRNA species###
##Record the length distribution of those walking distance###
##As a control, find the nearset A from the 3' end of piRNA species###

BEGIN { unshift @INC,"/home/ww74w/git/smallRNA_analysis/Utils/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;



open IN, "/home/ww74w/smallRNApipeline/pipeline_dm/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while (my $line=<IN>)
{
    if($line=~/>(.+) type/)
    {
        $chr="chr$1";
    }
    else
    {
        chomp $line;
        $genome{$chr}=$line;
    }
}
close(IN);

my %chrsize;
open IN, "/home/ww74w/smallRNApipeline/pipeline_dm/common/dm3.chrom.sizes";
while(my $line=<IN>) {
    next if($line=~/#/);
    chomp $line;
	@l=split(/\t/,$line);
    $chr=$l[0];
	$chrsize{$chr}=$l[1];
}
close(IN);

open IN, "$ARGV[0]";
open OUT, ">$ARGV[0].$ARGV[1].out";
 while(<IN>) #complementarity level
{
    chomp;
    split(/\t/);
    next if (length($_[4])>29 || length($_[4])<23);
    #$total{$file}+=$_[5]/$_[6];
    next if (/data/);
   
    for($i=1;$i<=$chrsize{$chr};$i++)
    {
        if($_[3] eq '+')
        {
            $start=$_[2]+$i-1;
            $str=substr($genome{$_[0]},$start,1);
            if($str eq $ARGV[1])
            {
            $step{$str}{$i}+=$_[5]/$_[6];
            last;
            }
        }
        else
        {
            $start=$_[1]-$i+1; #for bed format file, the end is the actural end coordinate+1 see UCSC bed format
            $str=substr($genome{$_[0]},$start,1);
            $str=&revfa($str);
            if($str eq $ARGV[1])
            {
            $step{$str}{$i}+=$_[5]/$_[6];
            last;
            }
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