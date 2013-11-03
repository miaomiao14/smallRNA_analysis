#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
# rule is p1 and p17-21 doesn't need to pair but p2-10 need, the rest 1mm
# simplifized version with prefix 16nt
# input as norm.bed ( no header)

####################################################################################################################################
#revised by WEI WANG
#11/14/2011
# 5'-5' distance is 10nt
# ask for at least 9nt complementarity
# allow 1 mm at the 1st position of the guide piRNA
####################################################################################################################################


open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while(<IN>) { if (/>(.+) type/) { $chr="chr$1";} else { chomp; $genome{$chr}=$_;}}

# make bowtie index
for ($i=0; $i<@ARGV; $i++) {
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
$hash2{$file}{substr($_[4],0,10)}+=$_[5]/$_[6];
#for ($n=1;$n<=20;$n++) {
#if ($_[3] eq '+') { $start=$_[1]+$n-17; $str=substr($genome{$_[0]},$start,16); $str=&revfa($str); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
#else {  $start=$_[2]-$n; $str=substr($genome{$_[0]},$start,16); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
#}
#}
########################################################################################################################################
if ($_[3] eq '+') { $start=$_[1]-1; $str=substr($genome{$_[0]},$start,10); $str=&revfa($str); $hash1{$file}{$str}+=$_[5]/$_[6];}
else {  $start=$_[2]-10; $str=substr($genome{$_[0]},$start,10); $hash1{$file}{$str}+=$_[5]/$_[6];}}
#########################################################################################################################################

if ($total{$file}>10) {
open OUT, ">$file.seq";
foreach (keys %{$hash2{$file}}) { print OUT "$_\t$hash2{$file}{$_}\n" if (length($_)==10);}
########################################################################################################################################
open OUT, ">$file.ref.fa";
foreach (keys %{$hash1{$file}}) { print OUT ">$_\t$hash1{$file}{$_}\n$_\n" if (length($_)==10);$k++;}
`bowtie-build $file.ref.fa $file`;
#########################################################################################################################################
 } }


#open OUT1,">pp_fraction";
#print OUT1 "pp";
#print OUT1 map {"\t$_"} @ARGV;
#print OUT1 "\n";

#print "pp";
#print map {"\t$_"} @ARGV;
#print "\n";

# bowtie mapping and score calculating
for ($i=0; $i<@ARGV; $i++) {
   print $ARGV[$i];
   for ($j=0; $j<=$i; $j++) {
   $file1=fileparse($ARGV[$i]); 
   @namefield=split(/\./,$file1);
   $name1=$namefield[2]."_".$namefield[1];
   $file1=$name1;
   $file2=fileparse($ARGV[$j]);
   @namefield=split(/\./,$file2);
   $name2=$namefield[2]."_".$namefield[1];
   $file2=$name2;
if ($total{$file1}<10 || $total{$file2}<10) { print "\t-10";} 
else {
$X=0; $Z=0;
%score=();
$count_N=0;





##************************************************************#

`bowtie $file1 -r -a -v 1 -p 8 $file2.seq --suppress 1,4,6,7 | grep + > $file1.$file2.bowtie.out`;
  %NTM=();
  open IN, "$file1.$file2.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
   
   if ($_[3]=~/(\d):/) { next if ($1!=9);} #0 -based index#only allow mismatch at the first position of the guide piRNA, corresponding to the 10th position of the target piRNA
   $NTM{$_[2]}++;
   }
  open IN, "$file1.$file2.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
  if ($_[3] eq "") {$t_9_nt=substr($_[2],9,1); $g_0_nt=&revfa($t_9_nt); $score{$g_0_nt.$t_9_nt}+=$hash1{$file1}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};}
  elsif ($_[3]=~/(\d):(\w)>(\w)/)
  { next if ($1!=9);  # no seed mismatches 0 based
  $g_0_nt=&revfa($2);
  $score{$g_0_nt.$3}+=$hash1{$file1}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};}}

##************************************************************#

open OUT2, ">$file1.$file2.frequency";
#print OUT2 "guide piRNA from $file2\n";
#print OUT2 "target piRNA from $file1\n";
  %NTM=();
 `bowtie $file2 -r -a -v 1 -p 8 $file1.seq --suppress 1,4,6,7 | grep + > $file2.$file1.bowtie.out`;
  open IN, "$file2.$file1.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
   
   if ($_[3]=~/(\d):/) { next if ($1!=9);}
   $NTM{$_[2]}++;}

   open IN, "$file2.$file1.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
   if ($_[3] eq "") {$t_9_nt=substr($_[2],9,1); $g_0_nt=&revfa($t_9_nt); $score{$g_0_nt.$t_9_nt}+=$hash1{$file2}{$_[1]}*$hash2{$file1}{$_[2]}/$NTM{$_[2]};}
   elsif ($_[3]=~/(\d):(\w)>(\w)/) { next if ($1!=9);
     $g_0_nt=&revfa($2);$score{$g_0_nt.$3}+=$hash1{$file2}{$_[1]}*$hash2{$file1}{$_[2]}/$NTM{$_[2]};}}
   foreach $Di_nt(keys %score)
   {

      print OUT2 "$file2\_$file1\t$Di_nt\t$score{$Di_nt}\n";
      print "$file2\_$file1\t$Di_nt\t$score{$Di_nt}\n";
      }
    foreach $Di_nt(keys %score)
   {
   $score{$Di_nt}=0;
   } 
     }
   

}
  #print "\n";
  #print OUT1 "\n";
  
}
 `rm *.ebwt`;
