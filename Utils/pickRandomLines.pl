#!/usr/bin/perl


##WEI WANG
##07/05/2012
## USAGE <INPUT> <COUNT lines>
use File::Basename;
use List::Util qw(shuffle);

$n=` wc -l $ARGV[0] | cut -d ' ' -f1`;
$n=$n-1;

open IN, "$ARGV[0]";
if(defined $ARGV[2])
{
	open OUT, ">$ARGV[2]/$ARGV[0].random.$ARGV[1].lines";
}
else
{
	open OUT, ">$ARGV[0].random.$ARGV[1].lines";
}

if($n<=100000000)
{
	my @lines = shuffle(<IN>);
	print OUT @lines[0 .. $ARGV[1]];
}
else
{
	my @lines = shuffle 0..$n ;
	$k=0;
	$i=0;
	$count=$ARGV[1];
	@line_num=@lines[0..$count-1];
	@line_num=sort { $a <=> $b } @line_num;
	$line_num=shift @line_num;
	while ($l=<IN>)
	{
		if($k<$count)
		{
			if ($line_num ==$i)
			{
			print OUT $l;
			$line_num=shift @line_num;
			$k++;
			}
		}
		$i++;
	}

}

close(OUT);
close(IN);
