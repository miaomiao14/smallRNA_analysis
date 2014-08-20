#!/usr/bin/perl
#03/03/2014
#WEI WANG
#check the shared 5 end between SRA and DEG

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;
use Getopt::Long;
use Cwd;

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "restrict_digts.pm";
require "parseGZInputFile.pm";
require "parseInputFile.pm";


my $numArgs = $#ARGV + 1;
my $format='normbed';
my $CLONE='SRA';
my $winSize=100;
my $inFileName1='';
my $inFileName2='';
my $outputDir=getcwd;
my $help="Help";
my $ok = GetOptions( 'help|?' => \$help,'inputFile1|i=s' => \$inFileName1,'inputFile2|j=s' => \$inFileName2, 'outputDir|o=s' => \$outputDir, 'format|F=s' => \$format, 'Clone|C=s' => $CLONE, 'windowSize|w=i' => \$winSize );

if($numArgs < 5 || !$ok || !$inFileName1 || !$outputDir )
{
  Usage();
  exit;
}


		my $file1="";	#$file[i] is the target strand,$file[j] is the guide strand
		$file1=fileparse($inFileName1);
		$file1 =~ /(.*)\.bed*.*/;
		my $filename1="";
		$filename1=$1;
    
		my $file2=""; #file2 is the guide strand
		$file2=fileparse($inFileName2);
		$file2 =~ /(.*)\.bed*.*/;
		my $filename2="";
		$filename2=$1;
   		
   		open PPSCORE, ">$outputDir/${filename1}_${filename2}.shared5end.stat.txt";
   		
		my %pp=(); my %pos=(); my %sharedpp=();my %sharedpos=();
					
		#read zip input file

		my %readMap1=();my %readMap2=();
		

		my $gz1 = gzopen($inFileName1, "rb") or die "Cannot open $inFileName1: $gzerrno\n" ;
		parseInputFile(\%readMap1,$gz1,$CLONE);
	    $gz1->gzclose();
			
		#read zip input file			
		my $gz2 = gzopen($inFileName2, "rb") or die "Cannot open $inFileName2: $gzerrno\n" ;
		parseInputFile(\%readMap2,$gz2,$CLONE);
	    $gz2->gzclose();
			
     	#calculate the PPscore for each pair of piRNAs with 5'-5' distance n	
		foreach my $strand (keys %readMap1)
		{
			if ($strand && exists $readMap2{$strand})
			{
				foreach my $chr (keys %{$readMap1{$strand}})
				{
					foreach my $cor (keys  %{$readMap1{$strand}{$chr}})
					{
						if( $readMap2{$strand}{$chr}{$cor})
						{
							$sharedpp{$cor}+=$pp{$cor};###todo
							$sharedpos{$cor}+=$pos{$cor};###todo
						}
					}
				}
			}
		}

		#pp
		my $totalppSpecies=0;
		my $totalppReads=0;
		$totalppSpecies=scalar (keys %pp);
		map { $totalppReads+=$_ } values %pp;
		
		my $sharedppSpecies=0;
		my $sharedppReads=0;
		$sharedppSpecies=scalar (keys %sharedpp);
		map { $sharedppReads+=$_ } values %sharedpp;
		
		$totalppReads=&restrict_num_decimal_digits($totalppReads,3);			
		$sharedppReads=&restrict_num_decimal_digits($sharedppReads,3);
		#fraction
		my $ppSpeciesF=0;
		$ppSpeciesF=$sharedppSpecies/$totalppSpecies;
		$ppSpeciesF=&restrict_num_decimal_digits($ppSpeciesF,3);
		my $ppReadsF=0;
		$ppReadsF=$sharedppReads/$totalppReads;
		$ppReadsF=&restrict_num_decimal_digits($ppReadsF,3);
		#pos
		my $totalposSpecies=0;
		my $totalposReads=0;
		$totalposSpecies=scalar (keys %pos);
		map { $totalposReads+=$_ } values %pos;
		
		my $sharedposSpecies=0;
		my $sharedposReads=0;
		$sharedposSpecies=scalar (keys %sharedpos);
		map { $sharedposReads+=$_ } values %sharedpos;
		#fraction
		my $posSpeciesF=0;
		$posSpeciesF=$sharedposSpecies/$totalposSpecies;
		$posSpeciesF=&restrict_num_decimal_digits($posSpeciesF,3);
		my $posReadsF=0;
		$posReadsF=$sharedposReads/$totalposReads;									
		$posReadsF=&restrict_num_decimal_digits($posReadsF,3);
		$totalposReads=&restrict_num_decimal_digits($totalposReads,3);			
		$sharedposReads=&restrict_num_decimal_digits($sharedposReads,3);
		
		print PPSCORE "$sharedppSpecies\t$totalppSpecies\t$ppSpeciesF\t$sharedppReads\t$totalppReads\t$ppReadsF\t$sharedposSpecies\t$totalposSpecies\t$posSpeciesF\t$sharedposReads\t$totalposReads\t$posReadsF\n";
		


close(PPSCORE);


sub usage
{
	print "\nUsage:$0\n\n\t";
	print "REQUIRED\n\t";
	print "-i  <inputfile1>\n\t";
	print "-j  <inputfile2 [optional]>\n\t";
	print "-o  <outputdir>\n\t";
	print "-f  <input format [bed|bedscore|normbed]>\n\t";
	#print "-t  <transposon name[optional]>\n\t";
	print "This perl script is to calculate the shared 5'end between smallRNAs and degradome reads\n\t";
	print "The input is bed file which includes both +(or sense) and -(or antisense) mappers\n\t";
	print "It's maintained by WEI WANG. If you have any questions, please contact wei.wang2\@umassmed.edu\n";
	
	exit(1);
	}

