#!/usr/bin/perl


##Find the nearest U from the 3' end of piRNA species###
##Record the length distribution of those walking distance###
##As a control, find the nearset A from the 3' end of piRNA species###

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;

###extract the sequences by bedtools already
###from bedtools nuc
###
my $inFile=$ARGV[0];
my $pattern=$ARGV[1];
my $patLen=length($pattern);
my $col=$ARGV[2]; #column # where the sequence is; 16
my %stepBySpe=();
my %stepByRed=();
open IN, "$inFile";

while(my $line=<IN>) #complementarity level
{
	next if($line=~/^#/);
    chomp $line;
    my @l=split(/\t/,$line);
 	my $seq=$l[$col];
	for ($i=0; $i<length($seq);$i++)
	{
		my $step=$i+1;
		my $charac=substr($seq,$i,$patLen);
		if(lc $charac eq lc $pattern)
		{
			$stepBySpe{$pattern}{$step}++;
			$stepByRed{$pattern}{$step}+=$l[3]/$l[4];
			last;
		}
		
	}
    
}
close(IN);

open OUT, ">$inFile.first.$pattern.steps.txt";
foreach my $str (keys %stepBySpe)
{
    foreach my $distance (sort {$a <=> $b } keys %{$stepBySpe{$str}} )
    {
        print OUT "$str\t$distance\t$stepBySpe{$str}{$distance}\t$stepByRed{$str}{$distance}\n";
    }
}
close(OUT);
