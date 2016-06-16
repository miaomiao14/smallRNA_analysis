#!/usr/bin/perl
## get read counts from bw file for the specified region
## modified from jie's code fetchBWScore.pl
## 2013-03-01, tushikui@gmail

use POSIX;   # for floor and ceil function
die "perl $0 <bigwig> <anchor bed> <tssUpLen> <nPtsGenebody> <ttsDnLen> <tss_tts_step> \n" if @ARGV<6;

my $bw_file     =$ARGV[0];   # input bigwig data file
my $genebed     =$ARGV[1];   # anchor bed (gene info)
my $tssUpLen    =$ARGV[2];   # define upstream length for TSS
my $nPtsGenebody=$ARGV[3];   # num of points from gene body
my $ttsDnLen    =$ARGV[4];   # define downstream length for TTS
my $tss_tts_step=$ARGV[5];   # interval length (resolution) for up/down-stream


open GIN, $genebed;
my( $extStart, $extEnd, $n1, $n2, $n3 );
while( my $t_line=<GIN> )
{	
	# each line: chr chrStart chrEnd genename val strand
	chomp $t_line;
	my @geneinfo = split /\t/, $t_line;

	$chr    = @geneinfo[0];   # chromosom
	$strand = @geneinfo[5];   # strand
	
	if ( $strand eq "+" ) {
		$extStart = @geneinfo[1] - $tssUpLen;
		$extEnd   = @geneinfo[2] + $ttsDnLen;
		$n1       = floor( $tssUpLen/$tss_tts_step );
		$n3       = floor( $ttsDnLen/$tss_tts_step );
	}
	if ( $strand eq "-" ) {
		$extStart = @geneinfo[1] - $ttsDnLen;
		$extEnd   = @geneinfo[2] + $tssUpLen;
		$n1       = floor( $ttsDnLen/$tss_tts_step );
		$n3       = floor( $tssUpLen/$tss_tts_step );
	}
	next if ( $extStart < 0 );

	# get reads pileup values
	my $x1 = `bigWigSummary $bw_file $chr $extStart @geneinfo[1] $n1 2>/dev/null`;
	my $x2 = `bigWigSummary $bw_file $chr $geneinfo[1] $geneinfo[2] $nPtsGenebody 2>/dev/null`;
	my $x3 = `bigWigSummary $bw_file $chr $geneinfo[2] $extEnd $n3 2>/dev/null`;
	chomp $x1;  chomp $x2;  chomp $x3;

	# check
	$x1 = &check_empty_and_strand( $x1,$strand,$n1 );
	$x2 = &check_empty_and_strand( $x2,$strand,$nPtsGenebody );
	$x3 = &check_empty_and_strand( $x3,$strand,$n3 );

	if ( $strand eq "+" ) { print $x1,"\t",$x2,"\t",$x3,"\n"; }
	if ( $strand eq "-" ) { print $x3,"\t",$x2,"\t",$x1,"\n"; }
}
close GIN;

sub check_empty_and_strand 
{
	my( $x,$strand,$n ) = @_;
	if ( $x eq "" ){ 
		join("\t",(0)x$n);
	}else{
		$x =~ s/n\/a/0/g;  # replace n/a with 0
		my @s = split( "\t",$x );
		@s = reverse @s if ( $strand eq "-" );
		join( "\t",@s );
	}
}

#my $bw=$ARGV[2];
#my $p=2*$bw/$ARGV[3];
#
#open anchor,$ARGV[1];
#while(<anchor>)
#{
#	chomp;
#	my ($chr,$pos,$strand)=(split)[0,1,5];
#	my $s=$pos-$bw;
#	my $e=$pos+$bw;
#	next if $s<0;
#	my $l=`bigWigSummary $ARGV[0] $chr $s $e $p 2>/dev/null`;
#	chomp $l;
#	my @s=();
#	if($l eq "")
#	{
#		@s=(0)x$p; 
#	}
#	else
#	{
#		$l=~s/n\/a/0/g;
#		@s=split("\t",$l);
#		@s=reverse @s if $strand eq "-";
#	}
#	print join("\t",@s),"\n";
#}
