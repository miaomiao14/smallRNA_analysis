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
use Compress::Zlib;
use Data::Dumper;
use threads;
use threads::shared;
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FileHandle;
use Benchmark qw(:hireswallclock);



my $numArgs = $#ARGV + 1;
#my $ok = getopts('i:o:F:D:w:', \%Options);

my $format='normbed';
my $CLONE='SRA';
my $winSize=100;

my $inFileName='';
my $outDir=getcwd;
my $help="Help";
my $ok = GetOptions( 'help|?' => \$help,'inputFile|i=s' => \$inFileName, 'outputDir|o=s' => \$outDir, 'format|F=s' => \$format, 'Clone|C=s' => $CLONE, 'windowSize|w=i' => \$winSize );

if($numArgs < 5 || !$ok || !$inFileName || !$outDir )
{
  Usage();
  exit;
}

my $file1=fileparse($inFileName);
my %readMap; #a hash to store strand, chr, coordinates(one end)

my $fh = FileHandle->new;
if ($fh->open("< $inFileName")) {
	parseInputFile(\%readMap,$fh);
    $fh->close;
}

open ARROUT, ">$outDir/temp.mt.array.sorted.txt";
print ARROUT Dumper(%readMap);
close(ARROUT);

###processing the distance 
my %disScoreChr=();
my %disScore=();
###for mutex
share(%disScoreChr);
#share(%disScore);
####variables for threads
my $nb_process = 10;
my @running = ();
my @Threads;

foreach my $strand (keys %readMap)
{
	if (!exists($disScoreChr{$strand})) {
		my %anon : shared;
		$disScoreChr{$strand} = \%anon;
	}
	foreach my $chr (keys %{$readMap{$strand}})
	{
		if (!exists($disScoreChr{$strand}{$chr})) {
			my %anon1 : shared;
			$disScoreChr{$strand}{$chr} = \%anon1;
		}
		&disProcess( \%{$readMap{$strand}{$chr}},$strand,$chr,\%disScoreChr);
        # my $trd= threads -> new ( sub {disProcess( \%{$readMap{$strand}{$chr}},$strand,$chr,\%disScoreChr) } );
        # push (@Threads, $trd);
        # my $tid = $trd->tid;
        # print "  - starting thread $tid\n";
	}#chr
}#strand

# foreach my $thr (@Threads)
# {
# 	print $thr->tid, " threadID\n";
# 	$thr->join;
# 	print "joining...\n";
# }

#sum up all distance frequencies for different strand, chr
while (my ($key, $value) = each %disScoreChr) {
	while (my ($key1, $value1) = each %{$value}) {
		while (my ($key2, $value2) = each %{$value1}) {
			$disScore{$key2}+=$value2;
		}
	}
}
# foreach my $strand (keys %disScoreChr)
# {
# 	foreach my $chr (keys %{$disScoreChr{$strand}})
# 	{
# 		foreach my $dis (keys %{$disScoreChr{$strand}{$chr}})
# 		{
# 			$disScore{$dis}+=$disScoreChr{$strand}{$chr}{$dis};
# 		}
# 	}
# }


#write the output to files

my $OUTDIR=$outDir;
# open OUT, ">$OUTDIR/$file1.5-5.distance.distribution" or die "cannot write to $OUTDIR: $!\n";
# foreach my $strand (keys %disScoreChr)
# {
# 	foreach my $chr (sort keys %{$disScoreChr{$strand}})
# 	{
# 		foreach my $dis (sort { $a <=> $b } keys %{$disScoreChr{$strand}{$chr}})
# 		{
# 			print OUT "$strand\t$chr\t$dis\t$disScoreChr{$strand}{$chr}{$dis}\n";
# 		}
# 	}
# }
# close (OUT);

open OUT, ">$OUTDIR/$file1.5-5.distance.distribution.summary" or die "cannot write to $OUTDIR: $!\n";
foreach my $dis (sort { $a <=> $b } keys %disScore)
{
	print OUT "$dis\t$disScore{$dis}\n";
}
close(OUT);


sub disProcess
{
	my ($readMapChrRef,$strand,$chr,$disScoreChrRef) =@_;
	my @sorted5end = &sort_hash_key( $readMapChrRef ); ##sort by the numerical value of the key, put it into an array
    foreach  (my $k=0;$k< $#sorted5end;$k++)
    {
        my $dis=0;
        foreach (my $j=$k+1;$j< $#sorted5end ;$j++)
        {
            $dis=$sorted5end[$j]-$sorted5end[$k]; ##why the $dis smaller than 0?
            if($dis>0 && $dis <=$winSize)
            {
				if (! exists $disScoreChrRef->{$strand}{$chr}{$dis}) {
					# my $val : shared;
					$disScoreChrRef->{$strand}{$chr}{$dis} = 0;
				}
                # print "3";
				# my $minDis : shared = ();
				my $minDis = &min($readMapChrRef->{$sorted5end[$j]},$readMapChrRef->{$sorted5end[$k]});
				# print "4";
				# {
					$disScoreChrRef->{$strand}{$chr}{$dis}+=$minDis;
# 				}
				# print "5";
                #$disScoreRef->{$dis}+=$minDis;
            }
            else
            {
				# print "6";
                $j=scalar( @sorted5end );
            }
        }
    }
}



sub parseInputFile
{
  		my ($readMapRef, $fileHandle) = @_;
		


        	if($format eq "normbed")
        	{
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
				}#while
        	}
			if($format eq "bed")
			{
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
					my $ntmreads;
					
					($chr,$bedstart,$bedend,$seq,$ntmreads,$strand)= split(/\t/,$line);
					$reads=$ntmreads;
					$ntm=1;
					$len=$bedend-$bedstart;
					
					if($CLONE eq "SRA") #length range to filter for piRNAs; for DEG, no filtering should be processed
	        		{
						next if ($len>29 || $len<23);
	        		}
	        		#($reads,$ntm,$dep)=split(/,/,$l[3]);
	        		if ($strand eq '+') #strand information is not included in the data, but in the file name
	        		{
		            	$readMapRef->{$strand}{$chr}{$bedstart}+=$reads/$ntm;
		            	#$plus_end{$l[0]}{$l[1]}=$_[2];
		        	}
		        	else
		        	{
		            	$readMapRef->{$strand}{$chr}{$bedend}+=$reads/$ntm;
		        	}
				}#while
        	}
        	if($format eq "bed2")#bo's definition
			{
				
				while(my $line = $fileHandle->getline())
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
        	}#while		
        			
		}#if
} #sub
		





sub Usage
{
	print "\n=======================================================================\n"; 
	print "\nUSAGE:$0\n\n -i <inputFile> -o <outputDir> -F <input-file type> <bed|normbed|bed2> -D <SRA | DEG> -w <window size of distance>\n\n"; 
	print " [-i <file>]\tinput file containing tags/reads in normbed/BED/bed2 format\n";
	print " [-o <file>]\toutput file into which directory will be stored (default: current working directory)\n";
	print " [-F <str>]\tinput-file type, currently bed,bed2(customized) and normbed are supported\n\t\t(default: normbed)\n";
	print " [-C <str>]\tcloning method (SRA or DEG, length matters, default:SRA)\n";
	print " [-w <int>]\tscanning window size downstread of each signal (default: 100)\n";
 	print "This perl script is count the frequency of 5'-5'end distances of smallRNAs(23-29) from the same strand\n";
	print  "\n=======================================================================\n\n";
	exit(1); 
} 
