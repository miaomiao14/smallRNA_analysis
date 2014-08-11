#!/usr/bin/perl

#Author: Wei Wang;
#email: weiwangumms@gmail.com
#the distribution of 5'-5' end distance of piRNAs from the same strand.

#the first trial of perl multithreading
#08/11/2014



BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "Statstics.pm";
require "Jia.pm";
BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "sort_hash_key.pm";
use File::Basename;
use Compress::Zlib;

use threads;
use threads::shared;
use strict;
use warnings;
#use Getopt::Long;
use Getopt::Std;
my %Options;

my $numArgs = $#ARGV + 1;
my $ok = getopts('i:o:F:D', \%Options);
if($numArgs < 4 || !$ok || !exists($Options{i}) || !exists($Options{o}) ) 
{ 
  Usage();
  exit;
}

	my $inFileName = $Options{i} ;
	my $outDir = $Options{o};
	my $format='normbed';
	if ( exists($Options{F}) )
  	{ 
		$format=$Options{F};
	}
	my $CLONE='SRA';
	if ( exists($Options{D}) )
  	{
	 	$CLONE=$Options{D};
  	}
  	my $winSize=100;
  	if ( exists($Options{w}) )
  	{
	 	$winSize=$Options{w};
  	}

my $file1=fileparse($inFileName);  

###main function; main hashes to store the data
my %plus=(); 
my %plus_end=(); 
my %minus=();
my $mins_end=();

my %readMap; #a hash to store strand, chr, coordinates(one end)
my $gz="";
$gz = gzopen($inFileName, "rb") or die "Cannot open $inFileName: $gzerrno\n" ;
parseInputFile(\%readMap,$gz);
$gz->gzclose();

###processing the distance 

my %disScoreChr=();
my %disScore=();

foreach my $strand (keys %readMap)
{  
	foreach my $chr (keys %{$readMap{$strand}})
	{
		my $plus_chr_ref = \%{$readMap{$strand}{$chr}};
		my @sorted5end = &sort_hash_key( $plus_chr_ref ); ##sort by the numerical value of the key, put it into an array
        &disProcess(\@sorted5end,\%readMap,$strand,$chr,\%disScore,\%disScoreChr);
        

	}#chr     
}#strand


#write the output to files

my $OUTDIR=$outDir;   
open OUT, ">$OUTDIR/$file1.5-5.distance.distribution" or die "cannot write to $OUTDIR: $!\n";
foreach my $strand (keys %disScoreChr)
{
	foreach my $chr (sort keys %{$disScoreChr{$strand}})
	{
		foreach my $dis (sort { $a <=> $b } keys %{$disScoreChr{$strand}{$chr}})
		{
			print OUT "$strand\t$chr\t$dis\t$disScoreChr{$strand}{$chr}{$dis}\n";
		}
	}
}
close (OUT);
       
open OUT, ">$OUTDIR/$file1.5-5.distance.distribution.summary" or die "cannot write to $OUTDIR: $!\n";
foreach my $dis (sort { $a <=> $b } keys %disScore)
{
	print OUT "$dis\t$disScore{$dis}\n";
}
close(OUT);


sub parseInputFile
{
  		my ($readMapRef, $fileHandle) = @_;
		
		while($fileHandle->gzreadline(my $line) > 0)
		{ 
			next if ($line=~/data/);
			chomp $line; 
			$line=~s/\s+/\t/g; #change all spaces to tab
			
			my $chr;
			my $bedstart;
			my $bedend;
			my $strand;
			my $seq;
			my $reads;
			my $ntm;
			my $len;

        	if($format eq "normbed")
        	{
        		($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$line);
				$len=$bedend-$bedstart+1;
				
				if($CLONE eq "SRA") #length range to filter for piRNAs; for DEG, no filtering should be processed
        		{
					next if ($len>29 || $len<23);
        		}
				if ($strand eq '+')
				{
	            	$readMapRef->{$strand}{$chr}{$bedstart}+=$reads/$ntm; 
	            	#$plus_end{$chr}{$bedend}=$l[2];
	        	}
	        	else
	        	{
	            	$readMapRef->{$strand}{$chr}{$bedend}+=$reads/$ntm;
	        	}
        	}
			if($format eq "bed")
			{
				
				my ($chr,$bedstart,$bedend,$seq,$ntmreads,$strand)= split(/\t/,$line);
				$reads=$ntmreads;
				$ntm=1;
				$len=$bedend-$bedstart;
				
				if($CLONE eq "SRA") #length range to filter for piRNAs; for DEG, no filtering should be processed
        		{
					next if ($len>29 || $len<23);
        		}
        		#($reads,$ntm,$dep)=split(/,/,$l[3]);
        		if($file1=~/plus/i) #strand information is not included in the data, but in the file name
        		{
	            	$readMapRef->{$strand}{$chr}{$bedstart}+=$reads/$ntm;
	            	#$plus_end{$l[0]}{$l[1]}=$_[2];
	        	}
	        	if($file1=~/minus/i)
	        	{
	            	$readMapRef->{$strand}{$chr}{$bedend}+=$reads/$ntm;
	        	}
        	}
        	if($format eq "bed2")#bo's definition
			{
				($chr,$bedstart,$bedend,$reads,$ntm,$strand,$seq)= split(/\t/,$line);
				$len=$bedend-$bedstart;

        		if($CLONE eq "SRA")
        		{
					next if ($len>29 || $len<23);
        		}
        		if($strand eq "+") #strand information is not included in the data, but in the file name
        		{
	            	$readMapRef->{$strand}{$chr}{$bedstart}+=$reads/$ntm; 
	            	#$plus_end{$l[0]}{$l[1]}=$l[2];
	        	}
	        	else
	        	{
	            	$readMapRef->{$strand}{$chr}{$bedend}+=$reads/$ntm;
	        	}
        	}		
        			
		}#while
} #sub
		

        sub disProcess
        {
        	my ($sorted5endRef,$readMapRef,$strand,$chr,$disScoreRef,$disScoreChrRef) =@_;
        	
	        foreach  (my $k=0;$k< scalar( @{$sorted5endRef} );$k++)
	        {
	            my $dis=0;
	    
	            foreach (my $j=$k+1;$j<scalar( @{$sorted5endRef} );$j++)
	            {
	                $dis=${$sorted5endRef}[$j]-${$sorted5endRef}[$k]; ##why the $dis smaller than 0?
	                if($dis>0 && $dis <=$winSize)
	                {
	                    
	                    my $minDis=&min($readMapRef->{$strand}{$chr}{${$sorted5endRef}[$j]},$readMapRef->{$strand}{$chr}{${$sorted5endRef}[$k]});
	                    $disScoreChrRef->{$strand}{$chr}{$dis}+=$minDis;
	                    $disScoreRef->{$dis}+=$minDis;
	                }
	                else
	                {
	                    $j=scalar( @{$sorted5endRef} );
	                }
	            }
	    
	        }
        }



sub Usage
{
	print "\n=======================================================================\n"; 
	print "\nUSAGE:$0\n\n -i <inputFile> -o <outputDir> -F <input-file type> [bed|normbed|bed2] -D [SRA | DEG]\n\n"; 
	print " [-i <file>]\tinput file containing tags/reads in normbed/BED/bed2 format\n";
	print " [-o <file>]\toutput file into which directory will be stored\n";
	print " [-F <int>]\tinput-file type, currently bed,bed2(customized) and normbed are supported\n\t\t(default: normbed)\n";
	print " [-D <int>]\tcloning method (SRA or DEG, length matters)\n";
	print " [-w <int>]\tscanning window size downstread of each signal (default: 100)\n";
 	print "This perl script is count the frequency of 5'-5'end distances of smallRNAs(23-29) from the same strand\n";
	print  "\n=======================================================================\n\n";
	exit(1); 
} 
