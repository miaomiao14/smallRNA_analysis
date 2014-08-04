#use v5.10.1;
 #   my %hash = (red    => 1, blue   => 2, green  => 3,
  #              orange => 4, yellow => 5, purple => 6#,
  #              black  => 7, grey   => 8, white  => 9);
  #  my @array = qw(red blue green);
  #  say "some array elements in hash keys" if  @array ~~  %hash;
  #  say "some array elements in hash keys" if \@array ~~ \%hash;
  #  say "red in array" if "red" ~~  @array;
  #  say "red in array" if "red" ~~ \@array;
  #  say "some keys end in e" if /e$/ ~~  %hash;
  #  say "some keys end in e" if /e$/ ~~ \%hash;

use PDL;
use PDL::Char;
#$PDL::SHARE = $PDL::SHARE; # keep stray warning quiet
use Data::Dumper;
use strict;
use warnings;

my $ref= 'ATTCCGGG' ;

my @query= ('ATTGCGGG' ,'ATACCGGC'  );

my $source = PDL::Char->new($ref);
for my $str (@query) {
	my $match = PDL::Char->new($str);

	# Is a PDL::Char
	my $diff = $match == $source;
	#$diff = ! $diff;	
	# Show diff
	#print join(' ', $diff->list), "\n";
	print Dumper $diff->list,"\n";
	#my $d;
	#$d = Data::Dumper->new([$match, $source]);
	#$d->Indent(3);
	#print $d->Dump;
	#print "list\t",$diff,"\n";
	my @bits=join('',$diff->list);
	#print @bits,"\tbits\n";
	#print $lenofDiff,"\n";
	#print $diff->atstr(0,1),"\n";
	#my @bits=split(/''/,"$diff->list");
	foreach my $b ( $diff->list)
	#foreach my $b=( 1 << @bit)
	{
		if($b eq "1" or $b == 1 )
		{
			my $bscale=$b*2;
			print "$b\t",$bscale,"\n";
		}
		else
		{
			print 0,"\t0\n";
		}
	}
	#print "\n";
}
# my $char = PDL::Char->new( [['abc', 'def', 'ghi'], ['jkl', 'mno', 'pqr']] );
# print $char->atstr(2,1),"\n";
#print $char->atstr(1,0),"\n";
my $bits=unpack ("b*", 2);
print "$bits\n";

my $strtest="ACCCGG";
my @strarr=split(//,$strtest);
print Dumper $strtest,"\n";
foreach my $r (@strarr)
{
	print $r,"\t";
}
print "\n";


my %test=("AT"=>1,"CG"=>2);
my @matchedPairs=("AT","TA","CG","GC");
my $sum=0;
foreach my $p ( @matchedPairs )
{
	if(!$test{$p})
	{
	$sum+=1;
	}
}
my $hashValue= values %test;
print $hashValue, "\n"; 
#print "$sum\n";
