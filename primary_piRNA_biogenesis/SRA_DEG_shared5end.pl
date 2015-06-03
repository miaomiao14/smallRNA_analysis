#!/usr/bin/perl
#06/03/2015
#WEI WANG
#calculate the fractions of shared 5′end between SRA and DEG

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;
use Getopt::Long;
use Cwd;

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/utils/";}
require "restrict_digts.pm";
require "parseGZInputFile.pm";
require "parseInputFile.pm";


my $numArgs = $#ARGV + 1;
my $format1='bed2';
my $format2='bed2';
my $CLONE1='SRA';
my $CLONE2='DEG';
my $winSize=100;
my $inFileName1='';
my $inFileName2=''; ##reference strand
my $outputDir=getcwd;
my $help="Help";
my $ok = GetOptions( 'help|?' => \$help,'inputFile1|i=s' => \$inFileName1,'inputFile2|j=s' => \$inFileName2, 'outputDir|o=s' => \$outputDir, 'format1|f=s' => \$format1,'format2|g=s' => \$format2, 'Clone1|c=s' => $CLONE1,'Clone2|d=s' => $CLONE2, 'windowSize|w=i' => \$winSize );

if($numArgs < 5 || !$ok || !$inFileName1 || !$outputDir )
{
  usage();
  exit;
}


		my $file1="";	#
		$file1=fileparse($inFileName1);
		$file1 =~ /(.*)\.bed*.*/;
		my $filename1="";
		$filename1=$1;
    
		my $file2=""; #file2 is the guide strand
		$file2=fileparse($inFileName2);
		$file2 =~ /(.*)\.bed*.*/;
		my $filename2="";
		$filename2=$1;
   		

   		
		my %sharedRef=();my %sharedQuery=();
					
		#read zip input file

		my %readMapQuery=();my %readMapRef=();
		
		my $posFlag=1; ##an parameter passing on to the subroutine to direct how data are stored in hash
		my $gz1 = gzopen($inFileName1, "rb") or die "Cannot open $inFileName1: $gzerrno\n" ;
		parseInputFile(\%readMapQuery,$gz1,$format1,$CLONE1,$posFlag);
	    $gz1->gzclose();
			
		#read zip input file			
		my $gz2 = gzopen($inFileName2, "rb") or die "Cannot open $inFileName2: $gzerrno\n" ;
		parseInputFile(\%readMapRef,$gz2,$format2,$CLONE2,$posFlag);
	    $gz2->gzclose();
			
     	#calculate the PPscore for each pair of piRNAs with 5'-5' distance n	
		foreach my $pos (keys %readMapQuery)
		{
			if ($pos && exists $readMapRef{$pos})
			{
				$sharedQuery{$pos}+=$readMapQuery{$pos};
				$sharedRef{$pos}+=$readMapRef{$pos};

			}
		}

		#ref
		my $totalRefSpecies=0;
		my $totalRefReads=0;
		my $sharedRefSpecies=0;
		my $sharedRefReads=0;

		$totalRefSpecies=scalar (keys %readMapRef);
		map { $totalRefReads+=$_ } values %readMapRef;

		$sharedRefSpecies=scalar (keys %sharedRef);
		map { $sharedRefReads+=$_ } values %sharedRef;
		
		$totalRefReads=&restrict_num_decimal_digits($totalRefReads,3);			
		$sharedRefReads=&restrict_num_decimal_digits($sharedRefReads,3);
		#fraction
		my $refSpeciesF=0;
		$refSpeciesF=$sharedRefSpecies/$totalRefSpecies;
		$refSpeciesF=&restrict_num_decimal_digits($refSpeciesF,3);
		my $refReadsF=0;
		$refReadsF=$sharedRefReads/$totalRefReads;
		$refReadsF=&restrict_num_decimal_digits($refReadsF,3);
		#query
		my $totalquerySpecies=0;
		my $totalqueryReads=0;
		$totalquerySpecies=scalar (keys %readMapQuery);
		map { $totalqueryReads+=$_ } values %readMapQuery;
		
		my $sharedQuerySpecies=0;
		my $sharedQueryReads=0;
		$sharedQuerySpecies=scalar (keys %sharedQuery);
		map { $sharedQueryReads+=$_ } values %sharedQuery;
		#fraction
		my $querySpeciesF=0;
		$querySpeciesF=$sharedQuerySpecies/$totalquerySpecies;
		$querySpeciesF=&restrict_num_decimal_digits($querySpeciesF,3);
		my $queryReadsF=0;
		$queryReadsF=$sharedQueryReads/$totalqueryReads;									
		$queryReadsF=&restrict_num_decimal_digits($queryReadsF,3);
		$totalqueryReads=&restrict_num_decimal_digits($totalqueryReads,3);			
		$sharedQueryReads=&restrict_num_decimal_digits($sharedQueryReads,3);
		
		open PPSCORE, ">$outputDir/${filename1}_${filename2}.shared5end.stat.txt";
		print PPSCORE "sharedRefSpecies\ttotalRefSpecies\trefSpeciesF\tsharedRefReads\ttotalRefReads\trefReadsF\tsharedQuerySpecies\ttotalquerySpecies\tquerySpeciesF\tsharedQueryReads\ttotalqueryReads\tqueryReadsF\n";
		print PPSCORE "$sharedRefSpecies\t$totalRefSpecies\t$refSpeciesF\t$sharedRefReads\t$totalRefReads\t$refReadsF\t$sharedQuerySpecies\t$totalquerySpecies\t$querySpeciesF\t$sharedQueryReads\t$totalqueryReads\t$queryReadsF\n";
		


close(PPSCORE);


sub usage
{
	print "\nUsage:$0\n\n\t";
	print "REQUIRED\n\t";
	print "-i  <inputfile1>\n\t";
	print "-j  <inputfile2 [optional]>\n\t";
	print "-o  <outputdir>\n\t";
	print "-f  <inputfile1 format [bed|bed2|normbed]>\n\t";
	print "-g  <inputfile2 format [bed|bed2|normbed]>\n\t";
	print "-c  <inputfile1 clone method [SRA|DEG|RSQ]>\n\t";
	print "-d  <inputfile2 clone method [SRA|DEG|RSQ]>\n\t";
	print "This perl script is to calculate the shared 5'end between smallRNAs and degradome reads\n\t";
	print "The input is bed file which includes both +(or sense) and -(or antisense) mappers\n\t";
	print "It's maintained by WEI WANG. If you have any questions, please contact wei.wang2\@umassmed.edu\n";
	
	exit(1);
	}

