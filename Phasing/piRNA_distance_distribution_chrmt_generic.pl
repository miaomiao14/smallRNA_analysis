#!/usr/bin/perl

#Author: Wei Wang;
#email: weiwangumms@gmail.com
#the distribution of 5'-5' end distance of piRNAs from the same strand.

#the first trial of perl multithreading
#08/11/2014

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "Statstics.pm";

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "sort_hash_key.pm";
use File::Basename;
# use Compress::Zlib;
#use Data::Dumper;
use threads;
use threads::shared;
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FileHandle;
# use Benchmark qw(:hireswallclock);
use Hash::Merge qw( merge );

# my $timespent;

my $numArgs = $#ARGV + 1;
my $format='normbed';
my $CLONE='SRA';
my $winSize=100;
my $inFileName1='';
my $inFileName2='';
my $end1='';
my $end2='';
my $outDir=getcwd;
my $help="Help";
my $ok = GetOptions( 'help|?' => \$help,'inputFile1|i=s' => \$inFileName1, 'end1|e1=s' => \$end1, 'inputFile2|j=s' => \$inFileName2,'end2|e2=s' => \$end2, 'outputDir|o=s' => \$outDir, 'format|F=s' => \$format, 'Clone|C=s' => $CLONE, 'windowSize|w=i' => \$winSize );

if($numArgs < 5 || !$ok || !$inFileName1 || !$outDir )
{
  Usage();
  exit;
}

my @Threadsinput=();

my $file1=fileparse($inFileName1);
my %readMap1; #a hash to store strand, chr, coordinates(one end)


###zlib is conflict with multithreading
###can not read gzip files

	 #my $trd= threads -> new ( sub {parseInputFile(\%readMap1,$inFileName1,$end1)} );	
	 #push (@Threadsinput, $trd);
	 parseInputFile(\%readMap1,$inFileName1,$end1);


if(!$inFileName2)
{
	$inFileName2=$inFileName1;
	$end2=$end1;
}

my $file2=fileparse($inFileName2);
my %readMap2; #a hash to store strand, chr, coordinates(one end)
if($end2!=$end1 or ($file2 ne $file1))
{

		#my $trd= threads -> new ( sub {parseInputFile(\%readMap2,$inFileName2,$end2)} );
		#push (@Threadsinput, $trd);
		parseInputFile(\%readMap2,$inFileName2,$end2);
}
else
{
	%readMap2=%readMap1;
}


# foreach my $thr (@Threadsinput)
# {
# 	print $thr->tid, " threadID\n";
# 	$thr->join;
# }

# open ARROUT, ">$outDir/temp.mt.array.sorted.txt";
# print ARROUT Dumper(%readMap);
# close(ARROUT);

###processing the distance 
my %disScoreChr=();
my %disScore=();

####variables for threads
my $nb_process = 10;
my @running = ();
my @Threads=();

foreach my $strand (keys %readMap1)
{

	foreach my $chr (keys %{$readMap1{$strand}})
	{
		
       # my $trd= threads -> new ( sub {disProcess( \%{$readMap{$strand}{$chr}},$strand,$chr,\%disScoreChr) } ); #each thread will write to the same hash, not a good idea
	   # my $trd= threads -> new ( sub {disProcess( \%{$readMap{$strand}{$chr}},$strand,$chr,\%{$disScoreChr{$chr}}) } ); #not try yet
		
		# print "\n\nNow trying without threading:\n\n";
		#my $starttime2 = Benchmark->new;
		#&disProcess( \%{$readMap{$strand}{$chr}},$strand,$chr,\%disScoreChr);		
		# my $finishtime2 = Benchmark->new;
		# $timespent = timediff($finishtime2,$starttime2);
		# print "\nDone!\n$chr at $strand Spent ". timestr($timespent),"\n";
		
		if( exists $readMap2{$strand}{$chr} )
		{
			my $trd= threads -> new ( sub {disProcess( \%{$readMap1{$strand}{$chr}},\%{$readMap2{$strand}{$chr}}, $strand,$chr) } ); #return value
					    push (@Threads, $trd);
					    my $tid = $trd->tid;

			# my $resultHashRef=disProcess( \%{$readMap1{$strand}{$chr}},\%{$readMap2{$strand}{$chr}}, $strand,$chr) ;
	# 		%disScoreChr=%{merge(\%disScoreChr,$resultHashRef)};
		}
		

		
        #print "  - starting thread $tid\n";
		
	}#chr
}#strand

foreach my $thr (@Threads)
{
	#print $thr->tid, " threadID\n";
	my $resultHashRef=$thr->join;

	%disScoreChr=%{merge(\%disScoreChr,$resultHashRef)};
	#print "joining...\n";
}


foreach my $strand (keys %disScoreChr)
{
	foreach my $chr (keys %{$disScoreChr{$strand}})
	{
		foreach my $dis (keys %{$disScoreChr{$strand}{$chr}})
		{
			foreach my $nuc (keys %{$disScoreChr{$strand}{$chr}{$dis}})
			{
				$disScore{$dis}{$nuc}+=$disScoreChr{$strand}{$chr}{$dis}{$nuc};
			}
		}
	}
}

#write the output to files; separate strand, chr
my $OUTDIR=$outDir;
open OUT, ">$OUTDIR/$file1.$end1.to.$file2.$end2.distance.distribution" or die "cannot write to $OUTDIR: $!\n";
foreach my $strand (keys %disScoreChr)
{
	foreach my $chr (sort keys %{$disScoreChr{$strand}})
	{
		foreach my $dis (sort { $a <=> $b } keys %{$disScoreChr{$strand}{$chr}})
		{
			foreach my $nuc (sort keys %{$disScoreChr{$strand}{$chr}{$dis}})
			{
				print OUT "$strand\t$chr\t$dis\t$nuc\t$disScoreChr{$strand}{$chr}{$dis}{$nuc}\n";
			}
		}
	}
}
close (OUT);
#write the output to files; accumulative distance score 
open OUT, ">$OUTDIR/$file1.$end1.to.$file2.$end2.distance.distribution.summary" or die "cannot write to $OUTDIR: $!\n";
foreach my $dis (sort { $a <=> $b } keys %disScore)
{
	foreach my $nuc (sort keys %{$disScore{$dis}})
	{
		print OUT "$dis\t$nuc\t$disScore{$dis}{$nuc}\n";
	}
}
close(OUT);


sub disProcess
{
	# my ($readMapChrRef,$strand,$chr,$disScoreChrRef) =@_;
	my ($readMapChrRef1,$readMapChrRef2,$strand,$chr) = @_;
	my %disScoreChrRef=();
	my @sortedEndFile1 = &sort_hash_key( $readMapChrRef1 ); ##sort by the numerical value of the key, put it into an array
	if($strand eq "+" or $strand eq "plus")
	{
	    foreach  (my $k=0;$k< $#sortedEndFile1;$k++)
	    {
	        my $dis=0;
			my $startPos=$sortedEndFile1[$k];
			
	        foreach (my $j=0;$j<=$winSize ;$j++) #$j is $dis
	        {
				$dis=$j;
	            my $newPos=$startPos+$dis;
	            if($readMapChrRef2->{$newPos})
	            {
					my $nuc=$readMapChrRef1->{$startPos}->[0].$readMapChrRef2->{$newPos}->[0];
					if (! exists $disScoreChrRef{$strand}{$chr}{$dis}) 
					{
						$disScoreChrRef{$strand}{$chr}{$dis}{$nuc} = 0;
					}
					
					my $minReads = &min($readMapChrRef1->{$startPos}->[1],$readMapChrRef2->{$newPos}->[1]);
					$disScoreChrRef{$strand}{$chr}{$dis}{$nuc}+=$minReads;
	            }

	        }#j
	    }#k
	}
	
	if($strand eq "-" or $strand eq "minus")
	{
	    foreach  (my $k=$#sortedEndFile1;$k>0;$k--)
	    {
	        my $dis=0;
			my $startPos=$sortedEndFile1[$k];
	        foreach (my $j=0;$j<=$winSize ;$j++)
	        {
				$dis=$j;
	            my $newPos=$startPos-$dis;
	            if($readMapChrRef2->{$newPos})
	            {
					my $nuc=$readMapChrRef1->{$startPos}->[0].$readMapChrRef2->{$newPos}->[0];
					if (! exists $disScoreChrRef{$strand}{$chr}{$dis}) 
					{
						$disScoreChrRef{$strand}{$chr}{$dis}{$nuc} = 0;
					}
					
					my $minReads = &min($readMapChrRef1->{$startPos}->[1],$readMapChrRef2->{$newPos}->[1]);
					$disScoreChrRef{$strand}{$chr}{$dis}{$nuc}+=$minReads;
	            }

	        }#j
	    }#k
	}	
	return (\%disScoreChrRef);
}



sub parseInputFile
{
  		my ($readMapRef, $file, $anchor) = @_;
		my $fileHandle = FileHandle->new;
		if ($fileHandle->open("< $file"))
		{

        	if($format eq "normbed")
        	{
        		while(my $line = $fileHandle->getline())
				{ 
					next if ($line=~/data/ or $line=~/\#/);
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
        		
	        		($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$line);
					$bedstart=$bedstart-1;
					
					
					&tagRecord($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm,$readMapRef,$anchor);
					
					
				}#while
        	}
			if($format eq "bed")
			{
				while(my $line = $fileHandle->getline())
				{ 
					next if ($line=~/data/ or $line=~/\#/);
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
					my $ntmreads;
					
					($chr,$bedstart,$bedend,$seq,$ntmreads,$strand)= split(/\t/,$line);
					$reads=$ntmreads;
					$ntm=1;
					
					
					&tagRecord($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm,$readMapRef,$anchor);
				}#while
        	}
        	if($format eq "bed2")#bo's definition
			{
				
				while(my $line = $fileHandle->getline())
				{ 
					next if ($line=~/data/ or $line=~/\#/);
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
					($chr,$bedstart,$bedend,$reads,$ntm,$strand,$seq)= split(/\t/,$line);
					
	
					&tagRecord($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm,$readMapRef,$anchor);
        	}#while		
        			
		}#if
	}#filehandle
	$fileHandle->close();
} #sub
		

sub tagRecord					 
{
	my ($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm,$readMapRef,$anchor) = @_ ;
	my $len=$bedend-$bedstart;
	if($CLONE eq "SRA") #length range to filter for piRNAs; for DEG, no filtering should be processed
	{
		return if ($len>29 or $len<23);
	}
	if($anchor==5)
	{
		my $nuc=substr($seq,0,1);
		
		if ($strand eq '+')
		{
        	$readMapRef->{$strand}{$chr}{$bedstart}->[1]+=$reads/$ntm; 
			$readMapRef->{$strand}{$chr}{$bedstart}->[0]=$nuc;
        	#$plus_end{$chr}{$bedend}=$l[2];
    	}
    	else
    	{
			my $fiveEnd=$bedend-1;
        	$readMapRef->{$strand}{$chr}{$fiveEnd}->[1]+=$reads/$ntm;
			$readMapRef->{$strand}{$chr}{$fiveEnd}->[0]=$nuc;
    	}
	}
	if($anchor==3)
	{
		my $nuc=substr($seq,-1);
		if ($strand eq '+')
		{
			#this 3' end is closed
			my $threeEnd=$bedend-1;
        	$readMapRef->{$strand}{$chr}{$threeEnd}->[1]+=$reads/$ntm; 
			$readMapRef->{$strand}{$chr}{$threeEnd}->[0]=$nuc;
        	#$plus_end{$chr}{$bedend}=$l[2];
    	}
    	else
    	{
			#if minus strand, 
			my $threeEnd=$bedstart;
        	$readMapRef->{$strand}{$chr}{$threeEnd}->[1]+=$reads/$ntm;
			$readMapRef->{$strand}{$chr}{$threeEnd}->[0]=$nuc;
    	}
	}
}#sub



sub Usage
{
	print "\n=======================================================================\n"; 
	print "\nUSAGE:$0\n\n -i <inputFile1> -j <inputFile2> -e1 <file1 End of Sequences> -e2 <file2 End of Sequences> -o <outputDir> -F <bed|normbed|bed2> -C <SRA|DEG> -w <windowSize of distance>\n\n"; 
	print " [-i <file>]\tinput file1 containing tags/reads in normbed/BED/bed2 format\n";
	print " [-j <file>]\tinput file2 containing tags/reads in normbed/BED/bed2 format\n";
	print " [-e1 <str>]\twhich end[5|3] of input file1\n";
	print " [-e2 <str>]\twhich end[5|3] of input file2\n";
	print " [-o <file>]\toutput file into which directory will be stored (default: current working directory)\n";
	print " [-F <str>]\tinput-file type, currently bed,bed2(customized) and normbed are supported\n\t\t(default: normbed)\n";
	print " [-C <str>]\tcloning method (SRA or DEG, length matters, default:SRA)\n";
	print " [-w <int>]\tscanning window size downstread of each signal (default: 100)\n";
 	print "This perl script is count the frequency of 5'-5'end distances of smallRNAs(23-29) from the same strand\n";
	print  "\n=======================================================================\n\n";
	exit(1); 
} 
