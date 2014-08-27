#!/usr/bin/perl

my $mergeFile=$ARGV[0];
my $sharedFile=$ARGV[1];
my $ppFile=$ARGV[2];

open M, "$mergeFile" or die "can not open file $mergeFile: $!";
open S, "$sharedFile" or die "can not open file $sharedFile: $!";
open P, "$ppFile" or die "can not open file $ppFile: $!";

my %sharedSpecies=();
while(my $line=<S>)
{
	chomp $line;
	my ($spe, $reads)=split(/\t/,$line);
	$sharedSpecies{$spe}=$reads;
}

my %ppSpecies=();
while(my $line=<P>)
{
	chomp $line;
	my ($spe, $g1t10)=split(/\t/,$line);
	push @{$ppSpecies{$spe}},$g1t10;
}

open O, ">$mergeFile.tag" or die "can not open file $mergeFile.tag to write: $!";
open OP, ">$mergeFile.ppseq.txt" or die "can not open file $mergeFile.ppseq.txt to write: $!";
while (my $line=<M>)
{
	if ($line =~ m/,/){
	chomp $line; 
	my ($g2x,$g1,$g2y,$seq,$r,$sum)=split(/\t/,$line);
	my @speciesf=split(/,/,$seq);
	my @tag=();
	for(my $i=0; $i <= $#speciesf; $i++)
	{
		if($sharedSpecies{$speciesf[$i]})
		{
			push @tag,"shared; ";
		}
		elsif($ppSpecies{$speciesf[$i]})
		{
			print OP "$speciesf[$i]\n";
			push @tag, "PP:";
			for(my $j=0; $j< @{$ppSpecies{$speciesf[$i]}}; $j++)
			{
				push @tag,"($ppSpecies{$speciesf[$i]}->[$j]),";
			}
			push @tag, "; ";
		}
		else
		{
			push @tag, "uniq no PP; ";
		}
	}
	print O $line,"\t",@tag,"\n";
	}
	
	
}
close(OP);
close(O);
close(P);
close(S);
close(M);