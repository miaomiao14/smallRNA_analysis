use strict;
use warnings;
use feature 'say';
my $A;
my $T;
my $G;
my $C;
my $fileIN;
my $fileOUT;

open $fileIN,  '<',"basecount.nfasta" or die "can't open file basecount.nfasta for reading";
open $fileOUT, '>','basecount.out' or die "can't open file basecount.out for writing";

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
  
  $total
}

close $fileIN;
close $fileOUT;