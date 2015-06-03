#!/usr/bin/perl
# use File::Basename;
# use threads;
# use threads::shared;
BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "Statstics.pm";

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "sort_hash_key.pm";
use File::Basename;
use Compress::Zlib;
use Data::Dumper;
use threads;
use threads::shared;
use strict;
use warnings;
use Getopt::Long;
use Cwd;


my $numArgs = $#ARGV + 1;
#my $ok = getopts('i:o:F:D:w:', \%Options);

my $format='normbed';
my $CLONE='SRA';
my $winSize=100;

my $inFileName='';
my $outDir=getcwd;
my $help="Help";
my $ok = GetOptions( 'help|?' => \$help,'inputFile|i=s' => \$inFileName, 'outputDir|o=s' => \$outDir, 'format|F=s' => \$format, 'Clone|C=s' => $CLONE, 'windowSize|w=i' => \$winSize );

if($numArgs < 5 || !$ok || !$inFileName || !$outDir )
{
  # Usage();
  exit;
}

my $file1=fileparse($inFileName);  

###main function; main hashes to store the data
my %plus=(); 
my %plus_end=(); 
my %minus=();
my $mins_end=();

my %readMap; #a hash to store strand, chr, coordinates(one end)
my $gz="";
$gz = gzopen($inFileName, "rb") or die "Cannot open $inFileName: $gzerrno\n" ;
parseInputFile(\%readMap,$gz);
$gz->gzclose();


# open ARROUT, ">$outDir/temp.mt.array.sorted.txt";
# print OUT Dumper(%readMap);
# close(OUT);

###processing the distance 

my %disScoreChr=();
my %disScore=();
###for mutex
# share(%disScoreChr);
#share(%disScore);
####variables for threads
my $nb_process = 10;
my @running = ();
my @Threads;

my @array;
share(@array);
my %hash;
share(%hash);
for (my $i=0; $i < 30; $i++) 
{
 $array[$i] = 0;
}
my @threadsarr;

for (my $i=0; $i<20; $i++)
{
	# add_hash_key(\@array,$i,$i*$i);
	my $trd= threads -> new (sub { add_hash_key(\%hash,$i,$i*$i) } );
	push (@threadsarr, $trd);
}

foreach my $thr (@threadsarr)
{
	print $thr->tid, " threadID\n";
	# $thr->join;
	$thr->join;
	print "joining...\n";
}

# for (my $i=0; $i<20; $i++)
# {
# 	print "$i=>$array[$i]\n";
# }

for my $key (keys %hash) {
	print "$key=>$hash{$key}\n";
}

sub add_hash_key
{
 my ($arr, $key, $value) = @_;
 # ${$arr}[$key] = $value;
 # print $key,"\t",${$arr}[$key],"\n";
 $arr->{$key} = $value;
 # print $key,
}
my $OUTDIR=$outDir;   
open OUT, ">$OUTDIR/$file1.5-5.distance.distribution" or die "cannot write to $OUTDIR: $!\n";
foreach my $strand (keys %disScoreChr)
{
	foreach my $chr (sort keys %{$disScoreChr{$strand}})
	{
		foreach my $dis (sort { $a <=> $b } keys %{$disScoreChr{$strand}{$chr}})
		{
			print OUT "$strand\t$chr\t$dis\t$disScoreChr{$strand}{$chr}{$dis}\n";
		}
	}
}
close (OUT);
       
open OUT, ">$OUTDIR/$file1.5-5.distance.distribution.summary" or die "cannot write to $OUTDIR: $!\n";
foreach my $dis (sort { $a <=> $b } keys %disScore)
{
	print OUT "$dis\t$disScore{$dis}\n";
}
close(OUT);
