#!/usr/bin/perl
# use mapper2 normalized by NTM & NTA
# normalization factor

use File::Basename;
if (@ARGV<1) {print "USAGE: <genelist.pl <mapper2> <norm factor>";}
$MAPPER2=$ARGV[0];
$filename = fileparse( $MAPPER2 );

$filename =~ /Phil.*.xkxh\.(.*)\.mapper2/;
$feature=$1;

open IN, $ARGV[0];
while(<IN>) {
chomp; split(/\t/);
$_[5]=~ /(.+)/;
$l=length($_[0]);


if (!exists $gene_piRNA{$1}) { $gene_piRNA{$1}=0;}
if (!exists $gene_piRNA_S{$1}) { $gene_piRNA_S{$1}=0;}
if (!exists $gene_piRNA_A{$1}) { $gene_piRNA_A{$1}=0;}
$gene_piRNA{$1}+=$_[1]/$_[6]/$_[7];
$gene_piRNA_S{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
$gene_piRNA_A{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');

}
close IN;

#open OUT,">$ARGV[0].list"; 
print "gene\t${feature}_piRNA\t${feature}_piRNA_sense\t${feature}_piRNA_antisense\t${feature}_piRNA_sense_fraction\n";

foreach (keys %gene_piRNA)
{

# normalization
if ($ARGV[1]){

$gene_piRNA{$_}=$gene_piRNA{$_}/$ARGV[1];
$gene_piRNA_S{$_}=$gene_piRNA_S{$_}/$ARGV[1];
$gene_piRNA_A{$_}=$gene_piRNA_A{$_}/$ARGV[1];
}

if ($gene_piRNA{$_}>0) {
$piRNA_sense_fraction=$gene_piRNA_S{$_}/$gene_piRNA{$_};
}
else {
$piRNA_sense_fraction='NA';
}

print "$_\t$gene_piRNA{$_}\t$gene_piRNA_S{$_}\t$gene_piRNA_A{$_}\t$piRNA_sense_fraction\n";
} 
