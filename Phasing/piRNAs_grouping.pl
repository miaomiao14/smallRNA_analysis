#!/usr/bin/perl

#01/24/2014
#the phasing for piRNAs only starts with U
#input format is normbed by default; if it's bed format, use bedtools nuc to extract the sequences first 

use File::Basename;
use Compress::Zlib;
if(scalar(@ARGV)<4)
{
        usage();
}
my $inputfile=$ARGV[0];
my $type=$ARGV[1];
my $outdir=$ARGV[2];
my $base=$ARGV[3];

#chrX	16503782	16503809	+	TAAAAACACTTCGTTGGATTTACGCGAG	1	1


#$., which is Perl's internal variable that keeps track of your current record number.
# $/ and $\ which are the input and output record separators 
#local $, = ','; print @arr;
$infilename=fileparse($inputfile);
$filename=fileparse($inputfile);
$filename=~s/.gz//;
open OUT, ">$outdir/$filename.$base" or die "cannot open $outdir/$filename.$base to write.$!";
open NOOUT, ">$outdir/$filename.non.$base" or die "cannot open $outdir/$filename.non.$base to write.$!";
	if($infilename=~/gz/)
	{
		my $gz="";
		$gz = gzopen($inputfile, "rb") or die "Cannot open $inputfile: $gzerrno\n" ;
		while($gz->gzreadline(my $line) > 0)
		{ 
			next if ($line=~/data/);
			chomp $line;
			$line= s/\s+/\t/g;
			my($chr,$start,$end,$strand,$seq,$reads,$ntm)=split(/\t/,$line);
			my $five="";

			$five=substr($seq,0,1);
			
			if($five eq $base)
			{
				print OUT $line,"\n";
				
			}
			else
			{
				print NOOUT $line,"\n";
			}

		}
		$gz->gzclose();
	}
close(NOOUT);
close(OUT);


sub usage
{
        print "\nUsage:$0\n\n\t";
        print "REQUIRED\n\t";
        print "inputfile type(bed|normbed) outdir base\n";
        print "This perl script is to split piRNAs(23-29) by their 5' end identity to different files\n";

        exit(1);
}