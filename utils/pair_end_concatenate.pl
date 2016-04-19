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
	my $skipline12=<fq1>;
	my $skipline21=<fq2>;
	my $skipline22=<fq2>;
	
	my @name1l=split(/\s/,$name1);
	my @name2l=split(/\s/,$name2);
	
	if($name1l[0] eq $name2l[0])
	{
		my $rcseq2=&revfa($seq2);
		print "$name1l[0]","\t","concatenated\n";
		print $seq1,"*",$rcseq2,"\n";
	}
	

}
close(fq2);
close(fq1);