#!/usr/bin/perl

#if you run it on hpcc

if(scalar(@ARGV)<3)
{
        usage();
}

$CLONE=$ARGV[3];
BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "Statstics.pm";
require "Jia.pm";
BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "sort_hash_key.pm";
use File::Basename;
use Compress::Zlib;
use Data::Dumper;
# the distribution of 5'-5' end distance of piRNAs from the same strand.

   $file1=fileparse($ARGV[0]);  
   #@namefield=split(/\./,$file1);
   #$name1=$namefield[0]."_".$namefield[2];
   #$file1=$name1;
   
   $format=$ARGV[1];
   
   %plus=(); %plus_end=(); %minus=();$mins_end=();
   
#	if($file1=~/gz/)
#	{
		my $gz="";
		$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
		while($gz->gzreadline($line) > 0)
		{ 
			chomp $line; $line=~s/\s+/\t/g;@l=split(/\t/,$line);  
			next if ($line=~/data/);
        
        	if($format eq "normbed")
        	{
				if($CLONE eq "SRA")
        		{
				next if (length($l[4])>29 || length($l[4])<23);
        		}
				if ($l[3] eq '+')
				{
	            	$plus{$l[0]}{$l[1]}+=$l[5]/$l[6];
	            	$plus_end{$l[0]}{$l[1]}=$l[2];
	        	}
	        	else
	        	{
	            	$minus{$l[0]}{$l[2]}+=$l[5]/$l[6];
	        	}
        	}
			if($format eq "bed")
			{

        		($reads,$ntm,$dep)=split(/,/,$l[3]);
        		if($file1=~/plus/i) #strand information is not included in the data, but in the file name
        		{
	            	$plus{$l[0]}{$l[1]}+=$reads/$ntm;
	            	$plus_end{$l[0]}{$l[1]}=$_[2];
	        	}
	        	if($file1=~/minus/i)
	        	{
	            	$minus{$l[0]}{$l[2]}+=$reads/$ntm;
	        	}
        	}
        	if($format eq "bed2")#bo's definition
			{
				
				my ($chr,$bedstart,$bedend,$reads,$ntm,$strand,$seq)= split(/\t/,$line);
				$len=$bedend-$bedstart;
				
        		if($CLONE eq "SRA")
        		{
					next if ($len>29 || $len<23);
        		}
        		if($strand eq "+") #strand information is not included in the data, but in the file name
        		{
	            	$plus{$chr}{$bedstart}+=$reads/$ntm;
	            	$plus_end{$l[0]}{$l[1]}=$l[2];
	        	}
	        	else
	        	{
	            	$minus{$chr}{$bedend}+=$reads/$ntm;
	        	}
        	}		
        			
		}
		$gz->gzclose();
#	}
  open OUT, ">$ARGV[2]/temp.hash.plus.txt";
print OUT Dumper(%plus);
close(OUT);
 open OUT, ">$ARGV[2]/temp.hash.minus.txt";
print OUT Dumper(%minus);
close(OUT);
    
    foreach $chr (keys %plus)
    {
        $plus_chr_ref = \%{$plus{$chr}};
        @plus_sort = &sort_hash_key( $plus_chr_ref ); ##sort by the numerical value of the key
        
        foreach  ($k=0;$k<$#plus_sort;$k++)
        {
            $dis=0;
    
            foreach ($j=$k+1;$j<$#plus_sort;$j++)
            {
                $dis=$plus_sort[$j]-$plus_sort[$k]; ##why the $dis smaller than 0?
                if($dis>0 && $dis <=100)
                {
                    $hash_dis{'plus'}{$chr}{$dis}+=&min($plus{$chr}{$plus_sort[$j]},$plus{$chr}{$plus_sort[$k]});
                    $h_dis{$dis}+=&min($plus{$chr}{$plus_sort[$j]},$plus{$chr}{$plus_sort[$k]});
                }
                else
                {
                    $j=$#plus_sort+1;
                }
            }
    
        }
        
    }
    foreach $chr (keys %minus)
        #%plus_lendis=&lendis_dist(%plus_sort);     ##calculate the distance distribution from the same strand, from the same chromosome
    {                                                #the adjacent piRNAs should be within 50nt of the current one
        if($minus{$chr})
        {
            
            @minus_sort=&sort_hash_key(\%{$minus{$chr}});
            #%minus_lendis=&lendis_dist(%minus_sort);
            foreach  ($k=0;$k<$#minus_sort;$k++)
            {
            $dis=0;
    
            foreach ($j=$k+1;$j<$#minus_sort;$j++)
            {
                $dis=$minus_sort[$j]-$minus_sort[$k];
                if($dis>0 && $dis <=100)
                {
                    $hash_dis{'minus'}{$chr}{$dis}+=&min($minus{$chr}{$minus_sort[$j]},$minus{$chr}{$minus_sort[$k]});
                    $h_dis{$dis}+=&min($minus{$chr}{$minus_sort[$j]},$minus{$chr}{$minus_sort[$k]});
                }
                else
                {
                    $j=$#minus_sort+1;
                }
            }
    
            }
        }
    }
 	$OUTDIR=$ARGV[2];   
    open OUT, ">$OUTDIR/$file1.5-5.distance.distribution";
    foreach $strand (keys %hash_dis)
    {
        foreach $chr (sort keys %{$hash_dis{$strand}})
        {
            foreach $dis (sort { $a <=> $b } keys %{$hash_dis{$strand}{$chr}})
            {
            print OUT "$strand\t$chr\t$dis\t$hash_dis{$strand}{$chr}{$dis}\n";
            }
        }
    }
    close (OUT);
    
    
    open OUT, ">$OUTDIR/$file1.5-5.distance.distribution.summary";
    foreach $dis (sort { $a <=> $b } keys %h_dis)
    {
        print OUT "$dis\t$h_dis{$dis}\n";
    }
    close(OUT);

sub usage
{
        print "\nUsage:$0\n\n\t";
        print "REQUIRED\n\t";
        print "inputfile type(bed|normbed) outdir clone\n";
        print "This perl script is count the frequency of 5'-5'end distances of smallRNAs(23-29) from the same strand\n";

        exit(1);
}
