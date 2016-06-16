#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;
BEGIN { unshift @INC,"/diag/home/netashawang/git/smallRNA_analysis/utils";}
require "Jia.pm";


my $fastq1=$ARGV[0];
my $filename1=fileparse($fastq1);

my $fastq2=$ARGV[1];
my $filename2=fileparse($fastq2);

if($fastq1 !~ /gz/)
{
	open fq1, "<$fastq1";
	open fq2, "<$fastq2";

	while (!eof(fq1) and !eof(fq2)) {
		my $name1=<fq1>;
		chomp $name1;
		my $name2=<fq2>;
		chomp $name2;
		my $seq1=<fq1>;
		chomp $seq1;
		my $seq2=<fq2>;
		chomp $seq2;
		my $skipline11=<fq1>;
		my $qc1=<fq1>;
		chomp $qc1;
		my $skipline21=<fq2>;
		my $qc2=<fq2>;
		chomp $qc2;
	
		my @name1l=split(/\s/,$name1);
		my @name2l=split(/\s/,$name2);
	
		if($name1l[0] eq $name2l[0])
		{
			my $rcseq2=&revfa($seq2);
			print "$name1l[0]","\t","concatenated\n";
			print $seq1,"*",$rcseq2,"\n";
			print "+\n";
			print "$qc1","*","$qc2","\n";
		}
	

	}
	close(fq2);
	close(fq1);
}
else
{

	my $fq1 = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
	my $fq2 = gzopen($ARGV[1], "rb") or die "Cannot open $ARGV[1]: $gzerrno\n" ;
	while($fq1->gzreadline(my $name1) > 0 and $fq2->gzreadline(my $name2) > 0)
	{
		$fq1->gzreadline(my $seq1);
		$fq2->gzreadline(my $seq2);
		$fq1->gzreadline(my $skipline1);
		$fq2->gzreadline(my $skipline2);
		$fq1->gzreadline(my $qc1);
		$fq2->gzreadline(my $qc2);
		
		chomp $name1;
		chomp $name2;
		
		chomp $seq1;
		chomp $seq2;
		
		chomp $qc1;
		chomp $qc2;
		
		my @name1l=split(/\s/,$name1);
		my @name2l=split(/\s/,$name2);
	
		if($name1l[0] eq $name2l[0])
		{
			my $rcseq2=&revfa($seq2);
			print "$name1l[0]","\t","concatenated\n";
			print $seq1,"*",$rcseq2,"\n";
			print "+\n";
			print "$qc1","*","$qc2","\n";
		}
	
		
		
	}
	$fq1->gzclose() ;
	$fq2->gzclose() ;
}