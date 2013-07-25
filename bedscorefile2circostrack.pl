#!/usr/bin/env perl
##This is a perl script to combine all the information from transposon mapper to circos data track
##input: DIR where the files are
##type mean or max or min
##binsize 
##prefix of file name
use File::Basename;

if(scalar(@ARGV)<7)
{
        usage();
}


$dir=shift @ARGV;
$normFacType=shift @ARGV;
$strand=shift @ARGV;
$metrictype=shift @ARGV;
$binsize=shift @ARGV;
#From normbedmapper2circos.pl
#open SOUT, ">$outdir/$filename\.$normFacType\.sense\.circos\.bed";
#from Bigwig
#\${j}.\${a[\$k]}.mean.\${BINSIZE}.txt
@files=`ls $dir/*.$normFacType.$strand.*.$metrictype.$binsize.*`;

#\${PIPELINE_DIRECTORY}/bedscorefile2circostrack.pl \$DIR mean \$bin \$file
$inputfilename=shift @ARGV;
($file,$dir)=fileparse($inputfilename);
$inputformat=shift @ARGV;
if($inputformat eq "mapper2")
{
$file =~ /(.*)\.mapper2/;
$filename=$1;
}
if($inputformat eq "normbed")
{
$file =~ /(.*)\.norm/;
$filename=$1;
}

open CHR, "/home/wangw1/pipeline/common/dm3.chrom.sizes";
while($read=<CHR>)
{
	chomp $read;
	@a=split(/\t/,$read);
	$chrsize{$a[0]}=$a[1];
}
close(CHR);

open OUT, ">$dir/$filename.$normFacType.$strand.$metrictype.$binsize.circos.txt";

foreach $f (@files)
{
	open SCORE, "$f";
	$sco=<SCORE>; #each chromosome is a line
	chomp $sco;
	@scores=split(/\t/,$sco);
	my ($filename,$dir)=fileparse($f);
	print "$filename\n";
	$filename=~/A.*.(chr.+)\.$type\.bin\.txt/;
	#$r[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
	$chrom=$1;
	$sizechr=$chrsize{$chrom};
	$nbins=$sizechr/$binsize;
		for($i=0;$i<($nbins-1);$i++)
		{
	  	$start=$binsize*$i;
	  	$end=$binsize*($i+1)-1;
	 		if($scores[$i] && $scores[$i] ne "n/a")
	  		{
	  		$binscore=$scores[$i];
	  		}
	  		else
	  		{
	  		$binscore=0;
	  		}
	  	print OUT "$chrom\t$start\t$end\t$binscore\n";
		}
		if($i!=0)
		{
		$start=$binsize*($i);
		$end=$sizechr;
			if($scores[$i] && $scores[$i] ne "n/a")
	  		{
	  		$binscore=$scores[$i];
	  		}
	  		else
	  		{
	  		$binscore=0;
			}
	  		print OUT "$chrom\t$start\t$end\t$binscore\n";
		}
		else
		{
		$start=0;
		$end=$sizechr;
			if($scores[$i] && $scores[$i] ne "n/a")
	  		{
	  		$binscore=$scores[$i];
	  		}
	  		else
	  		{
	  		$binscore=0;
	  		}
	  	print OUT "$chrom\t$start\t$end\t$binscore\n";   
		}
}
close(SCORE);
close(OUT);

sub usage
{
        print "\nUsage:$0\n\n\t";
        print "REQUIRED\n\t";
        print "inputdir normType[nnc|seqDep] strand[sense|antisense] metrics[max|mean] binsize filename fileformat[normbed|mapper2]\n";
        print "This perl script is to convert the inputfile to circos track!\n";

        exit(1);
}