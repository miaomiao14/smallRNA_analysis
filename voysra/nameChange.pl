#!/usr/bin/perl

my $nameM=$ARGV[0];
#my $sampleName=$ARGV[1];
my $fileList=`ls *.fastq.gz`;
open (IN, "<$nameM");
my %nameMap=();
while (my $line=<IN>)
{	
	chomp $line;
	my ($wx,$mir)=split(/\t/,$line);
	$nameMap{$wx}=$mir;
}
close(IN);
my @files=split(/\n/,$fileList);
#print $fileList;
foreach my $f (@files)
{
	foreach my $nm (keys %nameMap)
	{
		if ($f =~ /$nm/)
		{
			my $new=$f;
			$new=~s/$nm/$nameMap{$nm}/;
			`mv $f $new`;
		}
	}
}