#!/usr/bin/perl -w

### Script: illumina2fastq
### Author: ngs-interest.blogspot.com
### Last Update: 7 September 2011

use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX;

=head1 NAME

 illumina2fastq

=head1 SYNOPSIS

 illumina2fastq -[options]

 -h - help
 -i - input file (qseq or export file)
 -o - output file (default: write to out.fastq)

 example:
 illumina2fastq -i input_qseq.txt    ### input_qseq.txt is the "qseq" file and write output to out.fastq
 illumina2fastq -i input_qseq.txt -o outputfile.fastq ### input_qseq.txt is the "qseq" file and write output to outputfile.fastq
 

=head1 DESCRIPTION

 This script aims to create FASTQ file from the "qseq.txt" or "export.txt" file.

=cut

### option variables
my $help;
my $inFile;
my $outFile;
my $phredOffset;

### initialize option
Getopt::Long::Configure ('bundling');

if ( !GetOptions ('h'=>\$help, 'i=s'=>\$inFile, 'o=s'=>\$outFile) || !defined($inFile) ) {
 if ($help) {
  pod2usage(-verbose => 2);
 } else {
  pod2usage(1);
 }
}

### Initialize all variables
if (!defined($outFile)) {
 $outFile = "out.fastq";
}

if (!defined($phredOffset)) {
 $phredOffset = 0;
}

### Printing Information
print STDERR "Input File\t\t: $inFile\n";
print STDERR "Output File\t\t: $outFile\n";

print STDERR "\n";
print "Creating $outFile file...\n";

### File Initialization
open (IN, $inFile) || die "Cannot open $inFile for reading!\n";
open (OUT, ">$outFile") || die "Cannot open $outFile for writing!\n";

while (<IN>) {
 chomp;

 my @f = split(/\t/);
 $f[8] =~ s/\./N/g;

 my $id = "$f[0]-$f[1]-$f[2]-$f[3]-$f[4]-$f[5]";
 
 print OUT "\@$id\n$f[8]\n+$id\n$f[9]\n";
}
close (IN);
close (OUT);

print STDERR "Finish writing $outFile\n";
