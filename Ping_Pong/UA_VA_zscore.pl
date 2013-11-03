#!/usr/bin/perl

BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;

#@pairs=("AT","TA","GC","CG","AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");
#@VA=("AA","CA","GA");
#@HC=("AC","CC","TC");
#@DG=("AG","GG","TG");
#@BT=("CT","GT","TT");

%pairs=(
    VA => ["AA","CA","GA"],
    HC => ["AC","CC","TC"],
    DG => ["AG","GG","TG"],
    BU => ["CT","GT","TT"],
    UA => ["TA"],
    AU => ["AT"],
    CG => ["CG"],
    GC => ["GC"]
        );

open IN, "$ARGV[0]";

while ($line=<IN>)
{
    chomp $line;
    @l=split(/\t/,$line);
    $s_n{$l[0]}{$l[1]}=$l[3];
    $n_of_p{$l[0]}{$l[1]}=$l[2];
    #$s_b{$l[1]}{$l[0]}=$l[2];
}

close(IN);

 %rawscore=();
 %score=();
 %count_N=();
 
for ($i=1;$i<=20;$i++)
{
    if($s_n{$i})
    {
        foreach $pair ( keys %pairs)
        {
            foreach $p (@{$pairs{$pair}})
            {
                $s_n{$i}{$p}=0 if (!exists $s_n{$i}{$p});
                #$rawscore{$i}{$pair}+=$s_n{$i}{$p};
                $score{$pair}{$i}+=$s_n{$i}{$p};
                #if($i==10)
		#{
			$ncount{$pair}{$i}+=$n_of_p{$i}{$p};
		#} 
                
            }
            print "$i\t$pair\t$score{$pair}{$i}\n";
            $count_N{$pair}++ if ($score{$pair}{$i}>0);
        }
    }
}



open OUT, ">$ARGV[0].zscore.out";
@field=split(/\./,$ARGV[0]);


foreach $p (keys %pairs)
{
    $X{$p}=$score{$p}{10};
    delete $score{$p}{10};
    $NofPairs=$ncount{$p}{10};
    $std=&std(values %{$score{$p}});
    if ($std>0 && $count_N{$p}>=5) { $Z{$p}=($X{$p}-&mean(values %{$score{$p}}))/$std;} else {$Z{$p}=-10;}
    print OUT "$field[0]\-$field[1]\t$p\t$NofPairs\t$Z{$p}\n"; 
}
close(OUT);
