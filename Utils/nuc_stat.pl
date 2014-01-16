use strict;
use warnings;
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

open $fileIN,  '<',"$ARGV[0]" or die "can't open file basecount.nfasta for reading";
open $fileOUT, '>','$ARGV[1]' or die "can't open file basecount.out for writing";

while ( my $seq = <$fileIN> ) {

  next if $seq =~ /^>/;
  $seq =~ s/\n//g;
  say $seq;

  my @dna = split //, $seq;

  foreach my $element ( @dna ) {
    $A++ if $element =~ m/A/;
    $T++ if $element =~ m/T/;
    $G++ if $element =~ m/G/;
    $C++ if $element =~ m/C/;
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
  
  
}

close $fileIN;
close $fileOUT;