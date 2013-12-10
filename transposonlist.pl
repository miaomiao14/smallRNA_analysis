#!/usr/bin/perl
# use mapper2 normalized by NTM & NTA
# normalization factor
if (@ARGV<1) {print "USAGE: <transposonlist.pl <mapper2> <norm factor>";}

open IN, $ARGV[0];
while(<IN>) {
chomp; split(/\t/);
$_[5]=~ /(.+)\./;
$l=length($_[0]);

if (!exists $transposon{$1}) { $transposon{$1}=0;}
if (!exists $transposon_S{$1}) { $transposon_S{$1}=0;}
if (!exists $transposon_A{$1}) { $transposon_A{$1}=0;}
if (!exists $transposon_siRNA{$1}) { $transposon_siRNA{$1}=0;}
if (!exists $transposon_siRNA_S{$1}) { $transposon_siRNA_S{$1}=0;}
if (!exists $transposon_siRNA_A{$1}) { $transposon_siRNA_A{$1}=0;}
if (!exists $transposon_piRNA{$1}) { $transposon_piRNA{$1}=0;}
if (!exists $transposon_piRNA_S{$1}) { $transposon_piRNA_S{$1}=0;}
if (!exists $transposon_piRNA_A{$1}) { $transposon_piRNA_A{$1}=0;}

$transposon{$1}+=$_[1]/$_[6]/$_[7];
$transposon_S{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
$transposon_A{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');
if ($l==21) {
$transposon_siRNA{$1}+=$_[1]/$_[6]/$_[7];
$transposon_siRNA_S{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
$transposon_siRNA_A{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');
}
if ($l<=29 && $l>=23) {
$transposon_piRNA{$1}+=$_[1]/$_[6]/$_[7];
$transposon_piRNA_S{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
$transposon_piRNA_A{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');
}
}
close IN;

print "transposon\ttotal\ttotal_sense\ttotal_antisense\ttotal_sense_fraction\tsiRNA\tsiRNA_sense\tsiRNA_antisense\tsiRNA_sense_fraction\tpiRNA\tpiRNA_sense\tpiRNA_antisense\tpiRNA_sense_fraction\n";

foreach (keys %transposon)
{

# normalization
if ($ARGV[1]){
$transposon{$_}=$transposon{$_}/$ARGV[1];
$transposon_S{$_}=$transposon_S{$_}/$ARGV[1];
$transposon_A{$_}=$transposon_A{$_}/$ARGV[1];
$transposon_siRNA{$_}=$transposon_siRNA{$_}/$ARGV[1];
$transposon_siRNA_S{$_}=$transposon_siRNA_S{$_}/$ARGV[1];
$transposon_siRNA_A{$_}=$transposon_siRNA_A{$_}/$ARGV[1];
$transposon_piRNA{$_}=$transposon_piRNA{$_}/$ARGV[1];
$transposon_piRNA_S{$_}=$transposon_piRNA_S{$_}/$ARGV[1];
$transposon_piRNA_A{$_}=$transposon_piRNA_A{$_}/$ARGV[1];
}
if ($transposon_siRNA{$_}>0) {
$siRNA_sense_fraction=$transposon_siRNA_S{$_}/$transposon_siRNA{$_};
}
else {
$siRNA_sense_fraction='NA';
}
if ($transposon_piRNA{$_}>0) {
$piRNA_sense_fraction=$transposon_piRNA_S{$_}/$transposon_piRNA{$_};
}
else {
$piRNA_sense_fraction='NA';
}


print "$_\t$transposon{$_}\t$transposon_S{$_}\t$transposon_A{$_}\t",$transposon_S{$_}/$transposon{$_},"\t$transposon_siRNA{$_}\t$transposon_siRNA_S{$_}\t$transposon_siRNA_A{$_}\t$siRNA_sense_fraction\t$transposon_piRNA{$_}\t$transposon_piRNA_S{$_}\t$transposon_piRNA_A{$_}\t$piRNA_sense_fraction\n";
} 


