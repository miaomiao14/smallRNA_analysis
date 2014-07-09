#!/usr/bin/env perl

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "restrict_digts.pm";
require "Jia.pm";
require "supplementalPairing.pm";
use File::Basename;
use Compress::Zlib;
use PDL;  #it's needed for initialize, it has sum functio inside, clash with the one in Util
use PDL::Char;
$PDL::SHARE = $PDL::SHARE; 

my $inputFile=$ARGV[0];
%transPairSuppReadsTotal=();

#open IN, $inputFile or die:"failed to open $inputFile: $!"; 
my $gz = gzopen($inputFile, "rb") or die "Cannot open $inputFile: $gzerrno\n" ;
while($gz->gzreadline($line) > 0)
{
	@l=split(/\t/,$line);
	my $piGuideSuppSeq=$l[0];	
	my $piTargetIndexSuppSeq=$l[2]; #supplemental sequence for target strand
	my $diffstr='';
	if(length($piTargetIndexSuppSeq) == length($piGuideSuppSeq)) #check if the length of guideSuppSeq is the same as that of the $piTargetIndexSuppSeq
	{
		my $source = PDL::Char->new($piTargetIndexSuppSeq);
		my $match = PDL::Char->new($piGuideSuppSeq);
	      				     							
		my $diff = $match == $source; #bit comparison gives bit results: 0, not eqaul; 1 equal
		$diffstr=join('',$diff->list);#convert the list to a string
	}
	else
	{
		for(my $j=0;$j<$supLen;$j++)
		{
			$diffstr=$diffstr.'0'; #for default value
							
			}
	}
	
	$transPairSuppReadsTotal{$diffstr}+=$corPairReads;
}
$gz->gzclose();

my $basep=10;

my ($scaledSpeciesRef,$scaledReadsRef)=&SuppComBitSum(\%transPairSuppSpeciesTotal,\%transPairSuppReadsTotal);#for n=10
for(my $position=0; $position< @{$scaledReadsRef};$position++)
{
	my $supPos=$position+$basep+1;				
	print PPSEQSUPPVECTOR "transPair\tall\t$basep\t10\t$supPos\t$scaledSpeciesRef->[$position]\t$scaledReadsRef->[$position]\n";
}
	       					