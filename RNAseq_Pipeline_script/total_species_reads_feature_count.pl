#!/usr/bin/perl

#input file is mapper2

use File::Basename;
$MAPPER2=$ARGV[0];
open IN, "$MAPPER2";
while(<IN>) {
chomp; split(/\t/);
$species+=1/$_[6]/$_[7]; $reads+=$_[1]/$_[6]/$_[7];
$species_strand{$_[3]}+=1/$_[6]/$_[7];$reads_strand{$_[3]}+=$_[1]/$_[6]/$_[7];
}

$filename = fileparse( $MAPPER2 );
$filename =~ /Phil.*.xkxh\.(.*)\.mapper2/;
$feature=$1;
 print "$filename\t$feature\t$species\t$species_strand{sense}\t$species_strand{antisense}\t$reads\t$reads_strand{sense}\t$reads_strand{antisense}\n";