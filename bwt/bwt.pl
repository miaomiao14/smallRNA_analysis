#! /usr/bin/env perl
#this perl script is for implementing BWT (Burrows-Wheeler-Transform)
#contact: wei.wang2@umassmed.edu

my $str = "googol";
my $ending = "#";
my %rotations=();
my %K=(); #character count
sub permutation
{
	my ($str,$ending,$rotationRef)=@_;
	my $t=$str.$ending;
	$rotationRef->{0}=$t;
	print "$t\n";
	my @T=split(//,$t);
	for(my $i=0; $i<$#T;$i++ )
	{
		$K{$T[$i]}=0; #initialize %K, get the all unique characters;
		my $j=$i+1;
		my $tail=substr($t,0,$j);
		my $head=substr($t,$j);
		$rotationRef->{$j}="$head$tail";
		print "$head$tail\n";
	}
	
}
&permutation($str, $ending,\%rotations);

$K{$ending}=0; #initialize %K for the ending character

##sort the permutations lexicographically
#what sort (quicksort or other sorting?) is called here?
#

my @OCC=();
my $i=0;
foreach my $index ( sort { $rotations{$a} cmp $rotations{$b} } keys %rotations)
{
	# print $index,"\t",$rotations{$index},"\n";
	my $l=substr($rotations{$index},-1);
	my $f=substr($rotations{$index},0,1);
	$OCC[$i]=$K{$l};
	$K{$l}++;
	#print "$C[$i]\t$l\n";
	print $i,"\t",$K{$l},"\n";
	$i++;
}
print "\n\n";
foreach my $character (sort {$a cmp $b } keys %K)
{
	print $character,"\t","$K{$character}","\n";
}
print "\n\n";
foreach my $k (@OCC)
{
	print $k,"\t",$OCC[$k],"\n";
}
	

