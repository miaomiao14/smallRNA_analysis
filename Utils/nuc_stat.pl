#!/usr/bin/perl

#01/16/2014
#wei.wang2@umassmed.edu

#count the frequencies of nucleotide from sequences from clusters and transposons

#the input file is bed format
#use the bedtools nuc extract the sequences and Profiles the nucleotide content of intervals in a fasta file

#the input file is the output of bedtools nuc files with seq on the last column.

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;
use feature 'say';


my $A;
my $T;
my $G;
my $C;
my $total;
my $A_fre;
my $T_fre;
my $G_fre;
my $C_fre;
my $fileIN;
my $fileOUT;

open $fileIN,  "$ARGV[0]" or die "can't open file basecount.nfasta for reading";
open $fileOUT, ">$ARGV[1]" or die "can't open file basecount.out for writing";
open OUT, ">$ARGV[2]" or die "can't open file basecount.out for writing";

while ( my $line = <$fileIN> ) {

  next if $line =~ /^>/;
  chomp $line;

  my($chr,$start,$end,$name,$temp,$strand,$ATperc,$GCperc,$A1,$C1,$G1,$T1,$N,$O,$Len,$seq0)=split(/\t/,$line);
  my $seq = uc $seq0;
  my $rbound=$end-23; #assume the length of piRNAs are at least 23 nt long
  my $char = 'T';
  my $offset = 0;
  my $tindex = index($seq, $char, $offset);

  while ($tindex != -1) {

	

    my $newstart=$start+$tindex;
    if($newstart>$rbound)
    {
    	last;
    }
    else
    {
    my $string=substr($seq,$tindex,23);
    my $newend=$newstart+23;
    $newstart=$newstart+1;#to accomondate to norm.bed format
	print OUT "$chr\t$newstart\t$newend\t\+\t$string\t1\t1\n";
    $offset = $tindex + 1;
    $tindex = index($seq, $char, $offset);
    }

  }
  
  my $lbound=$start+22;
  $char = 'A';
  $offset = 0;
  $tindex = index($seq, $char, $offset);

  while ($tindex != -1) {

	

    my $newend=$start+$tindex;
    if($newend<$lbound)
    {
    	last;
    }
    else
    {
    my $string=substr($seq,$tindex,23);
    my $stringrc=&RevComp($string);
    my $newstart=$newend-23;
    $newstart=$newstart+1;#to accomondate to norm.bed format
	print OUT "$chr\t$newstart\t$newend\t\-\t$stringrc\t1\t1\n";
    $offset = $tindex + 1;
    $tindex = index($seq, $char, $offset);
    }

  }
  
  #say $seq;

  my @dna = split //, $seq;

  foreach my $element ( @dna ) {
    $A++ if $element =~ m/A/;
    $T++ if $element =~ m/T/;
    $G++ if $element =~ m/G/;
    $C++ if $element =~ m/C/;
  }
  

  
}

  say $fileOUT "A=$A";
  say $fileOUT "T=$T";
  say $fileOUT "G=$G";
  say $fileOUT "C=$C";
  
  $total=$A+$T+$C+$G;
  $A_fre=$A/$total;
  $C_fre=$C/$total;
  $G_fre=$G/$total;
  $T_fre=$T/$total;
  
  
  say $fileOUT "A_fre=$A_fre";
  say $fileOUT "T_fre=$T_fre";
  say $fileOUT "G_fre=$G_fre";
  say $fileOUT "C_fre=$C_fre";

close $fileIN;
close $fileOUT;
close(OUT);

my %RevCompBasePairs=qw/A T T A G C C G a t t a g c c g U A u a R Y r y Y R y r M K m k K M k m S S s s W W w w H D D H h d d h B V V B b v v b N N n n/;

sub RevComp {
  my $s=shift @_;
  my $new_s='';
  for my $b (split//,$s) {
   if(!exists $RevCompBasePairs{$b}) { $new_s.="?"; }
   else { $new_s=$RevCompBasePairs{$b}.$new_s; }
  }
  return $new_s;
}