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
    #XA => ["AA","CA","GA","UA"],
    VA => ["AA","CA","GA"],
    HC => ["AC","CC","TC"],
    DG => ["AG","GG","TG"],
    BU => ["CT","GT","TT"],
    #UX => ["TA","TC","TG","TT"],	
    UA => ["TA"],
    #UB => ["TC","TG","TT"],
    AU => ["AT"],
    #AX => ["AT","AC","AG","AA"],
    #XU => ["AT","CT","GT","TT"],
    CG => ["CG"],
    GC => ["GC"]
        );

open IN, "$ARGV[0]";

%s_n=();
%n_of_p=();
$total_s=0;
$total_r=0;
while ($line=<IN>)
{
    chomp $line;
    @l=split(/\t/,$line);
    $s_n{$l[0]}{$l[1]}=$l[3];    #number of reads
    $total_r+=$l[3];
    $n_of_p{$l[0]}{$l[1]}=$l[2]; #number of species
    $total_s+=$l[2];
    #$s_b{$l[1]}{$l[0]}=$l[2];
}

close(IN);

 %rawnReads=();
 %nReads=();
 %nSpecies=();
 open OUT, ">$ARGV[1]/$ARGV[0].raw.pair.nReads.out";
foreach $i (keys %s_n) 
{        foreach $pair ( keys %pairs)
        {
            foreach $p (@{$pairs{$pair}})
            {
                $s_n{$i}{$p}=0 if (!exists $s_n{$i}{$p});
                #$rawnReads{$i}{$pair}+=$s_n{$i}{$p};
                $nReads{$i}{$pair}+=$s_n{$i}{$p};
				$nSpecies{$i}{$pair}+=$n_of_p{$i}{$p};
	
                
            }
            print OUT "$i\t$pair\t$nSpecies{$i}{$pair}\t$nReads{$i}{$pair}\n";
            print "$i\t$pair\t$nSpecies{$i}{$pair}\t$nReads{$i}{$pair}\n";
            #$count_N{$pair}++ if ($nReads{$pair}{$i}>0);
        }

}
close(OUT);

#
#open OUT, ">$ARGV[1]/$ARGV[0].nReads.out";
##@field=split(/\./,$ARGV[0]);
#
#print OUT "PP pairs\tUA_reads_probability\tUA_pairs_probability\tAU_reads_probability\tAU_pairs_probability\n";
#foreach $i (keys %nReads)
#
#{
#    #$reads{UA}=$nReads{$i}{UA}*$total_r/($nReads{$i}{XA}*$nReads{$i}{UX});
#    #$species{UA}=$nSpecies{$i}{UA}*$total_s/($nSpecies{$i}{XA}*$nSpecies{$i}{UX});
#    $reads{UA}=$nReads{$i}{UA}/($nReads{$i}{UX});
#    $species{UA}=$nSpecies{$i}{UA}/($nSpecies{$i}{UX});
#    #$reads{AU}=$nReads{$i}{AU}*$total_r/($nReads{$i}{XU}*$nReads{$i}{AX});
#    #$species{AU}=$nSpecies{$i}{AU}*$total_s/($nSpecies{$i}{XU}*$nSpecies{$i}{AX});
#    $reads{AU}=$nReads{$i}{AU}/($nReads{$i}{AX});
#    $species{AU}=$nSpecies{$i}{AU}/($nSpecies{$i}{AX});
#    print OUT "$i\t$reads{UA}\t$species{UA}\t$reads{AU}\t$species{AU}\n"; 
#}
#close(OUT);
