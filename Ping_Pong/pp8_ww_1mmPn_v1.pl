#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
# rule is p1 and p17-21 doesn't need to pair but p2-10 need, the rest 1mm
# simplifized version with prefix 16nt
# input as norm.bed ( no header)


open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while(<IN>) { if (/>(.+) type/) { $chr="chr$1";} else { chomp; $genome{$chr}=$_;}}
$m=$ARGV[3];
# make bowtie index
for ($i=0; $i<$ARGV[2]; $i++) {
my $file=fileparse($ARGV[$i]);
    @namefield=split(/\./,$file);
    $name=$namefield[2]."_".$namefield[1];
    push @argos, $name;
    $file=$name;
open IN, $ARGV[$i];
while(<IN>) { chomp; split(/\t/);  
next if (length($_[4])>29 || length($_[4])<23);
next if (/data/);
$total{$file}+=$_[5]/$_[6];
$hash2{$file}{substr($_[4],0,16)}+=$_[5]/$_[6];
for ($n=1;$n<=20;$n++) {
if ($_[3] eq '+') { $start=$_[1]+$n-17; $str=substr($genome{$_[0]},$start,16); $str=&revfa($str); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
else {  $start=$_[2]-$n; $str=substr($genome{$_[0]},$start,16); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
}
}
if ($total{$file}>10) {
open OUT, ">$file.seq";
foreach (keys %{$hash2{$file}}) { print OUT "$_\t$hash2{$file}{$_}\n" if (length($_)==16);}
for ($n=1;$n<=20;$n++) {
open OUT, ">$file.ref.$n.fa";
foreach (keys %{$hash1{$file}{$n}}) { print OUT ">$_\t$hash1{$file}{$n}{$_}\n$_\n" if (length($_)==16);$k++;}
`bowtie-build $file.ref.$n.fa $file.$n`;
#`rm $file.ref.$n.fa`;
 } }
}

open OUT1,">pp_fraction";
print OUT1 "pp";
print OUT1 map {"\t$_"} @argos;
print OUT1 "\n";
open O, ">zscore.out";
# bowtie mapping and score calculating
for ($i=0; $i<$ARGV[2]; $i++) {
   for ($j=0; $j<=$i; $j++) {
      
   $file1=fileparse($ARGV[$i]); 
   @namefield=split(/\./,$file1);
   $name1=$namefield[2]."_".$namefield[1];
   $file1=$name1;
   $file2=fileparse($ARGV[$j]);
   @namefield=split(/\./,$file2);
   $name2=$namefield[2]."_".$namefield[1];
   $file2=$name2;
   
   print O "$file1-$file2";
   print "$file1-$file2";
   print OUT1 "$file1-$file2";
   
if ($total{$file1}<10 || $total{$file2}<10) { print "\t-10"; print O "\t-10"; print OUT1 "\tNA";} 
else {

 $X=0; $Z=0; %score=();$count_N=0;
open OUT, ">$file1.$file2.pp";
open OUT2, ">$file1.$file2.ppseq";
foreach ($n=1;$n<=20;$n++) {
# file1 as ref
 `bowtie $file1.$n -r -a -v 1 -p 8 $file2.seq --suppress 1,4,6,7 | grep + > $file1.$file2.$n.bowtie.out`;
#  foreach ($m=0;$m<=9;$m++)
#  {
      %NTM=();
      open IN, "$file1.$file2.$n.bowtie.out";
      while(<IN>) { chomp; split(/\t/);  
      if ($_[3]=~/(\d+):/) {
 
      next if ($1!=$m);}  # no seed mismatches 0 based
      $NTM{$_[2]}++; ##if there is 0 mm, there will be addition here
      #$count_N++;
      }
  
      open IN, "$file1.$file2.$n.bowtie.out";
      while(<IN>) { chomp; split(/\t/);
      if ($_[3]=~/(\d+):/) {
         next if ($1!=$m);}  # no seed mismatches 0 based
         $score{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};
         print OUT2 "NA\t$_[2]\n" if ($n==10);
      }
#  }
# `rm $file1.$file2.$n.bowtie.out`;  
# file2 as ref
  %NTM=undef;
 `bowtie $file2.$n -r -a -v 1 -p 8 $file1.seq --suppress 1,4,6,7 | grep + > $file2.$file1.$n.bowtie.out`;
 open IN, "$file2.$file1.$n.bowtie.out";
  while(<IN>) { chomp; split(/\t/);  
  if ($_[3]=~/(\d+):/) {
   next if ($1!=$m);}  # no seed mismatches 0 based
   $NTM{$_[2]}++;
  # $count_N++;
  }
   open IN, "$file2.$file1.$n.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
  if ($_[3]=~/(\d+):/) {
      next if ($1!=$m);}  # no seed mismatches 0 based
      $score{$n}+=$hash1{$file2}{$n}{$_[1]}*$hash2{$file1}{$_[2]}/$NTM{$_[2]};
      print OUT2 "$_[2]\tNA\n" if ($n==10);
  }
# `rm $file2.$file1.$n.bowtie.out`;
     $score{$n}=0 if (!exists $score{$n});
     print OUT "$n\t$score{$n}\n";
     $count_N++ if ($score{$n}>0);
}
   $X=$score{10}; delete $score{10};
   $std=&std(values %score); 
   if ($std>0 && $count_N>=5) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
   print O "\t$m\t$Z";
   print "\t$m\t$Z";
 
   if ($Z!=-10) { print OUT1 "\t",$X/$total{$file1}/$total{$file2}*1000000;} else { print OUT1 "\tNA";}
   $N1=`match.pl $file1.$file2.ppseq $file1.seq | sumcol+ 2`; chomp($N1);
   `cut -f2 $file1.$file2.ppseq > $file1.$file2.ppseq.temp`;
   $N2=`match.pl $file1.$file2.ppseq.temp $file2.seq | sumcol+ 2`; chomp($N2);
   open OUT3, ">$file1.$file2.pp_reads_percentage";
   print OUT3 "$N1\t$total{$file1}\t",$N1/$total{$file1},"\t","$N2\t$total{$file2}\t",$N2/$total{$file2},"\n";
}
  print O "\n";
  print "\n";
  print OUT1 "\n";
}

}
`rm *.ebwt`;
