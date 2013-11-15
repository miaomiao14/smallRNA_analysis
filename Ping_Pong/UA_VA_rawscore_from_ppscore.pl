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
    ($piwipair,$overlap,$basepair,$spe,$reads)=split(/\t/,$line);
    $n_of_r{$overlap}{$piwipair}{$basepair}=$reads;    #number of reads
    $total_r+=$reads;
    $n_of_s{$overlap}{$piwipair}{$basepair}=$spe; #number of species
    $total_s+=$reads;
    #$s_b{$l[1]}{$l[0]}=$l[2];
}

close(IN);

 %nReads=();
 %nSpecies=();
open OUT, ">$ARGV[1]/$ARGV[0].raw.pair.nReads.out";
foreach $overlap (keys %n_of_s) 
{        
	foreach $piwipp (keys %{$n_of_s{$overlap}})
	{
		foreach $pair ( keys %pairs)
        {
            foreach $p (@{$pairs{$pair}})
            {
                $n_of_s{$overlap}{$piwipp}{$p}=0 if (!exists $n_of_s{$overlap}{$piwipp}{$p});
                #$rawnReads{$i}{$pair}+=$s_n{$i}{$p};
                $nReads{$overlap}{$piwipp}{$pair}+=$n_of_s{$overlap}{$piwipp}{$p};
				$nSpecies{$overlap}{$piwipp}{$pair}+=$n_of_r{$overlap}{$piwipp}{$p};
	
                
            }
            print OUT "$overlap\t$piwipp\t$pair\t$nSpecies{$overlap}{$piwipp}{$pair}\t$nReads{$overlap}{$piwipp}{$pair}\n"; 
        }
        
        #print "$i\t$pair\t$nSpecies{$i}{$pair}\t$nReads{$i}{$pair}\n";

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
