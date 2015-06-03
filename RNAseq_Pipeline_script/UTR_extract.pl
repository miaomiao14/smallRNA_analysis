#!/usr/bin/perl

#

use File::Basename;
use List::Util qw(sum);
use warnings;
use strict;

my $longpep=$ARGV[0];
my $bed=$ARGV[1];
my $outdir=$ARGV[2];
my $pepType=$ARGV[3];
if(scalar(@ARGV)<4)
{
        usage();
}


`[ ! -f $longpep.complete ] && grep "complete" $longpep > $longpep.complete `;
open PEP, "$longpep.complete" or die "can not find file $longpep.complete: $!";
my %cds;
while(my $line=<PEP>)
{
	chomp $line;
	$line =~ s/>//;
	$line =~ s/\s+/\t/g;
	my @fields = split ( /\t/, $line );
	my @chrinfo=();
	if ($pepType eq "long")
	{
		@chrinfo=split(/:/,$fields[4]);
	}
	if ($pepType eq "best")
	{
		@chrinfo=split(/:/,$fields[8]);
	}
	$chrinfo[1]=~ /(\d+)-(\d+)\(\+\)/;
	my $cdsStart=$1-1; #convert to 0-based
	my $cdsEnd=$2;
	push @{$cds{$chrinfo[0]}}, [$cdsStart,$cdsEnd]; #one transcript can have multiple protein coding gene candicates
	#print "$chrinfo[0]\t$cdsstart\t$cdsend\n";
}
close(PEP);

#convert gtf (from cufflinks output) to bed first, use bedClip and bedtools sort already!
open BED, "$bed" or die "can not find file $bed: $!";
my %transcript;
while (my $line=<BED>)
{
	chomp $line;
	my ($chr, $start, $end, $name, $score, $strand, $x0, $x1, $x2, $blocknum, $blocksize, $blockstart) = split ( /\t/, $line );
	@{$transcript{$name}} = [$chr, $start, $end, $strand, $blocknum, $blocksize, $blockstart] ;
}
close(BED);


my $filename = fileparse( $bed );
$filename =~ /^(.*)\.bed$/;
my $prefix=$1;
$prefix.="_";

my $pepfile = fileparse($longpep);
$pepfile =~ /^(.*)\.pep/;
$prefix.=$1;

open FUTR, ">$outdir/$prefix.5UTR.bed" or die "can not find file $outdir/$prefix.5UTR.bed to write: $!";
#open FUTRL,">$outdir/$prefix.5UTR.length.list" or die "can not find file $outdir/$prefix.5UTR.length.list: $!";

open TEST, ">$outdir/test.txt";
foreach my $record ( sort keys %cds )
{
	print TEST $record,"\t";

		my $transcriptName=$record;
		my $exonTotalLen = 0;
		print TEST "$transcriptName\n";
		#print TEST "$record\t$name\t$transcript{$name}[0]->[0]\t$transcript{$name}[0]->[1]\n"; 
		if ($transcript{$transcriptName})
		{
			my ($chr, $transcriptStart, $transcriptEnd, $strand, $blocknum, $blocksize, $blockstart)=@{$transcript{$transcriptName}[0]};#each transcript one record
			my @exonSizes = split(/,/,$blocksize);
			$exonTotalLen = sum (@exonSizes);
			my @exonStarts = split(/,/,$blockstart);
			my $sizeOfArray = scalar @exonStarts ;
			#print TEST "$chr, $transcriptStart, $transcriptEnd, $strand, $blocknum, $blocksize, $blockstart","\n";
			foreach my $cdsIsoform (@{$cds{$record}})
			{
				my($cdsStart, $cdsEnd)=($cdsIsoform->[0],$cdsIsoform->[1]);
				#print TEST "$cdsStart, $cdsEnd","\n";
				my $cdsLen = $cdsEnd - $cdsStart;
				
				if ($strand eq "+")
				{
					my $FutrStart = $transcriptStart; #0-based
					#my $FutrEnd = $transcriptStart + $cdsStart -1; #need to be converted into genomic coordinates
					my $FutrEnd = $FutrStart; #initialize the End
					my $FutrLength = $cdsStart;
					
					my $lenLeft = $FutrLength;
					
					my $i = 0;
					my @FutrSize=();
					my @FutrStartPos=();
					if ($lenLeft !=0) #some protein coding candidates do not have 5UTR
					{
						while($lenLeft>0)
						{
							$lenLeft = $lenLeft - $exonSizes[$i];
							push @FutrSize, $exonSizes[$i];
							push @FutrStartPos, $exonStarts[$i];
							$i++;
						}
						my $j = $i;
						
						$i = $i - 1;
						$lenLeft = $lenLeft + $exonSizes[$i];
						pop @FutrSize;
						push @FutrSize, $lenLeft;
						$FutrEnd = $FutrEnd+$exonStarts[$i]+$lenLeft;
						print FUTR "$chr\t$FutrStart\t$FutrEnd\t$transcriptName\t0\t$strand\t$FutrEnd\t$FutrEnd\t0\t$j\t",join(",",@FutrSize),",\t",join(",",@FutrStartPos),",\n";
					}
				}
				else
				{
					my $FutrEnd = $transcriptEnd;
					my $FutrStart = $exonStarts[$sizeOfArray-1]; #the end of the first exon in the antisense gene
					
					my $FutrLength = $cdsStart;
					my $lenLeft = $exonTotalLen - $FutrLength;
					
					my $i = 0 ;
					my @FutrSize = @exonSizes;
					my @FutrStartPos= @exonStarts;
					
					if ($lenLeft != $exonTotalLen && $lenLeft !=0 ) #some protein coding candidates do not have 5UTR
					{
						my @FutrStartAdjustedPos=();
						while($lenLeft>0)
						{
							$lenLeft = $lenLeft - $exonSizes[$i];
							shift @FutrSize;
							shift @FutrStartPos;
							$i++;	
						}
						
						
						
						my $j =  scalar @FutrSize;
						if($lenLeft<0)
						{
							my $FutrSize0= -$lenLeft ; # if $lenLeft<0
							unshift @FutrSize, $FutrSize0;
							$j = scalar @FutrSize;
							
							$i = $i - 1;
							$lenLeft = $lenLeft + $exonSizes[$i];
							$FutrStart = $transcriptStart + $exonStarts[$i] + $lenLeft;
							@FutrStartAdjustedPos = map {$_-($exonStarts[$i] + $lenLeft)} @FutrStartPos; #adjust the exon start sites
							
							unshift @FutrStartAdjustedPos, 0;
						}
						else # if $lenLeft=0
						{
							#$i = $i - 1;
							#$lenLeft = $lenLeft + $exonSizes[$i];
							$FutrStart = $transcriptStart + $exonStarts[$i] ;
							@FutrStartAdjustedPos = map {$_-($exonStarts[$i])} @FutrStartPos;
							#shift @FutrStartAdjustedPos;
							#unshift @FutrStartAdjustedPos, 0;
						}
						print FUTR "$chr\t$FutrStart\t$FutrEnd\t$transcriptName\t0\t$strand\t$FutrEnd\t$FutrEnd\t0\t$j\t",join(",",@FutrSize),",\t",join(",",@FutrStartAdjustedPos),",\n";
						
					}
					
				}
				
			}
		}
}
close(TEST);
close(FUTR);
#close(FUTRL);

sub usage
{
        print "\nUsage:$0\n\n\t";
        print "REQUIRED\n\t";
        print "PEP(transdecoder longest_orf.pep) BED(cufflink_gtf_to_bed) OUTDIR PEPTYPE(long or best)\n";
        print "This perl script is to extract 5UTR of transcripts based on the annotation of transdecoder\n";

        exit(1);
}
