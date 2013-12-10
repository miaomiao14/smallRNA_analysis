#!/usr/bin/perl
# use mapper2 normalized by NTM & NTA
# normalization factor

use Compress::Zlib;
if (@ARGV<1) {print "USAGE: <transposonlist.pl <mapper2> <norm factor>";}

#open IN, $ARGV[0];
#while(<IN>) {

$gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[0]: $gzerrno\n" ;
while($gz->gzreadline($_) > 0)
{
	
chomp; split(/\t/);
$_[5]=~ /(.+)\./;
$l=length($_[0]);

if (!exists $transposon{$1}) { $transposon{$1}=0;}
if (!exists $transposon_S{$1}) { $transposon_S{$1}=0;}
if (!exists $transposon_A{$1}) { $transposon_A{$1}=0;}


$transposon{$1}+=$_[1]/$_[6]/$_[7];
$transposon_S{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'sense');
$transposon_A{$1}+=$_[1]/$_[6]/$_[7] if ($_[3] eq 'antisense');

}
#close IN;
$gz->gzclose();
print "transposon\tpiRNA\tpiRNA_sense\tpiRNA_antisense\tpiRNA_sense_fraction\n";
foreach (keys %transposon)
{

# normalization
if ($ARGV[1]){
$transposon{$_}=$transposon{$_}/$ARGV[1];
$transposon_S{$_}=$transposon_S{$_}/$ARGV[1];
$transposon_A{$_}=$transposon_A{$_}/$ARGV[1];

}


print "$_\t$transposon{$_}\t$transposon_S{$_}\t$transposon_A{$_}\t",$transposon_S{$_}/$transposon{$_},"\n";
}

