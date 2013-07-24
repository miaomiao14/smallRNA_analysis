#!/usr/bin/env perl


use File::Basename;
  
if(scalar(@ARGV)<1)
{
        usage();
}


if($ARGV[1])
{
	open CHR, "$ARGV[1]";
}
else
{
	open CHR, "/home/wangw1/ppeline/common/dm3.chrom.sizes";
}
while($read=<CHR>)
{
	chomp $read;
	@a=split(/\t/,$read);
	$chrsize{$a[0]}=$a[1];
}
if($ARGV[1])
{
	open NF, "$ARGV[2]";
}
else
{
	open NF, "/home/wangw1/pipeline/common/nf_2012_01.txt";
}
while($l=<NF>)
{
	chomp $l;
	@l=split(/\t/,$l);
	$nf{$l[0]}=$l[1];
}

open IN, "$ARGV[0]";
open PLUSOUT, ">$ARGV[0]\.plus\.circos.bed";
open MINUSOUT, ">$ARGV[0]\.minus\.circos.bed";
$file=fileparse($ARGV[0]);
$filename=`basename $file .xkxh.norm.bed`;
$normfactor=$nf{$filename}/1000000;
while($r=<IN>) 
{ 
	chomp $r;
	@r= split(/\t/,$r);
	next if ($r[5]>$chrsize{$r[2]});
	$start=$r[1]-1;
	if($r[4] eq "+")
	{
		print PLUSOUT "$r[0]\t$start\t$r[2]\t$r[5],$r[6],$normfactor\n";
	}
	else
	{
		print MINUSOUT "$r[0]\t$start\t$r[2]\t$r[5],$r[6],$normfactor\n";
	}
}

close(PLUSOUT);
close(MINUSOUT);
close(IN);

sub usage
{
        print "\nUsage:$0\n\n\t";
        print "REQUIRED\n\t";
        print "inputfile(norm.bed) CHRSIZE NormFactorFile\n";
        print "This perl script is convert the normbed file to bed file with 4th column is the score(reads,NTM,nf)\n";

        exit(1);
}