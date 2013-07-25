#!/usr/bin/env perl


use File::Basename;
use Compress::Zlib; 
if(scalar(@ARGV)<1)
{
        usage();
}

my $inFile=shift @ARGV;
my $type=shift @ARGV;
################################################################
my $chrSizeFile=shift @ARGV;
if($chrSizeFile)
{
	open CHR, "$chrSizeFile";
}
else
{
	open CHR, "/home/wangw1/ppeline/common/dm3.chrom.sizes";
}
while(my $chr=<CHR>)
{
	chomp $chr;
	my @chrLine=split(/\t/,$chr);
	$chrsize{$chrLine[0]}=$chrLine[1];
}
################################################################
my $normFacFile=shift @ARGV;
if($normFacFile)
{
	open NF, "$normFacFile";
}
else
{
	open NF, "/home/wangw1/pipeline/common/nf_2012_01.txt";
}
################################################################
while(my $l=<NF>)
{
	chomp $l;
	my @l=split(/\t/,$l);
	$nf{$l[0]}=$l[1];
}
################################################################
my $outdir=shift @ARGV;

$file=fileparse($inFile);
open IN, "$inFile";
$gz = gzopen($inFile, "rb") or die "Cannot open $inFile: $gzerrno\n" ;


if($type eq "normbed")
{
	$filename=`basename $file .norm.bed.gz`;
	$normfactor=$nf{$filename}/1000000;

	open PLUSOUT, ">$outdir/$file\.plus\.circos\.bed";
	open MINUSOUT, ">$outdir/$file\.minus\.circos\.bed";
	
	while($gz->gzreadline($r) > 0) #read zipped files
	#while($r=<IN>) 
	{ 
		chomp $r;
		@r= split(/\t/,$r);
		if (length($r[0])>=23 & length($r[0])<=34)
		{
			next if ($r[5]>$chrsize{$r[2]});
			$start=$r[1]-1;
			if($r[4] eq "+")
			{
				print PLUSOUT "$r[0]\t$start\t$r[2]\t$r[5],$r[6],$normfactor\n";
			}
			if($r[4] eq "-")
			{
				print MINUSOUT "$r[0]\t$start\t$r[2]\t$r[5],$r[6],$normfactor\n";
			}
		}
	}
	close(PLUSOUT);
	close(MINUSOUT);
	close(IN);
}
if($type eq "mapper2")
{

	$filename=`basename $file .mapper2.gz`;
	$normfactor=$nf{$filename}/1000000;
	open SOUT, ">$outdir/$file\.sense\.circos\.bed";
	open AOUT, ">$outdir/$file\.antisense\.circos\.bed";

	while($gz->gzreadline($r) > 0)	
	#while($r=<IN>) 
	{ 
		chomp $r;
		@r= split(/\t/,$r);
		if (length($r[0])>=23 & length($r[0])<=34)
		{
			$r[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
			next if ($3>$chrsize{$1});
			$ntm_nta=$r[6]*$r[7];
			$key="$1\t$2\t$3\t$r[1],$ntm_nta,$normfactor\n";
			if($r[3] eq "sense")
			{
				$hash{$key}=0;
			}
			else 
			{
				$hash{$key}=1;
			}
		}
	}

	foreach $line (keys %hash) 
	{
		if($hash{$line}==0) 
		{
			print SOUT "$line";
		}
		if($hash{$line}==1)
		{
			print AOUT "$line";
		} 
	}
	close(AOUT);
	close(SOUT);
	close(IN);
}
################################################################

sub usage
{
        print "\nUsage:$0\n\n\t";
        print "REQUIRED\n\t";
        print "inputfile(*.norm.bed|*.mapper2) type(normbed|mapper2) CHRSIZE NormFactorFile outdir\n";
        print "This perl script is convert the normbed file to bed file with 4th column is the score(reads,NTM,nf)\n";

        exit(1);
}