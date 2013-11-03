#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;


####################################################################################################################################
#revised by WEI WANG
#07/25/2012
# 5'-5' distance is 10nt
# ask for at least 15nt complementarity
# allow for 1 mm at the 10th position of the target piRNA
# Just record raw score
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
$hash2{$file}{substr($_[4],0,16)}+=$_[5]/$_[6];
#for ($n=1;$n<=20;$n++) {
#if ($_[3] eq '+') { $start=$_[1]+$n-17; $str=substr($genome{$_[0]},$start,16); $str=&revfa($str); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
#else {  $start=$_[2]-$n; $str=substr($genome{$_[0]},$start,16); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
#}
#}
########################################################################################################################################
if ($_[3] eq '+') { $start=$_[1]-7; $str=substr($genome{$_[0]},$start,16); $str=&revfa($str); $hash1{$file}{$str}+=$_[5]/$_[6];}
else {  $start=$_[2]-10; $str=substr($genome{$_[0]},$start,16); $hash1{$file}{$str}+=$_[5]/$_[6];}}
#########################################################################################################################################

if ($total{$file}>10) {
open OUT, ">$file.seq";
foreach (keys %{$hash2{$file}}) { print OUT "$_\t$hash2{$file}{$_}\n" if (length($_)==16);}
########################################################################################################################################
open OUT, ">$file.ref.fa";
foreach (keys %{$hash1{$file}}) { print OUT ">$_\t$hash1{$file}{$_}\n$_\n" if (length($_)==16);$k++;}
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

%score=();
%species=();

open OUT2, ">$file2.$file1.1.frequency";
#print OUT2 "$file1_$file2\t";



##************************************************************#

`bowtie $file1 -r -a -v 1 -p 8 $file2.seq --suppress 1,4,6,7 | grep + > $file2.$file1.bowtie.out`;
  %NTM=();
  open IN, "$file2.$file1.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
   
   if ($_[3]=~/(\d):/) { next if ($1!=0);} #0 -based index#only allow mismatch at the first position of the guide piRNA, corresponding to the 10th position of the target piRNA
   $NTM{$_[2]}++;
   }
  open IN, "$file2.$file1.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
  if ($_[3] eq "") {$g_0_nt=substr($_[2],0,1); $t_9_nt=&revfa($g_0_nt); $score{$g_0_nt.$t_9_nt}+=$hash1{$file1}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};
                    $species{$g_0_nt.$t_9_nt}{$_[2]}=1;}
  elsif ($_[3]=~/(\d):(\w)>(\w)/)
  { next if ($1!=0);  # no seed mismatches 0 based
  $t_9_nt=&revfa($2);
  $score{$3.$t_9_nt}+=$hash1{$file1}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};
  $species{$3.$t_9_nt}{$_[2]}=1;
  }
  }
foreach $Di_nt(keys %score)
{
   $n_of_species=scalar (keys %{$species{$Di_nt}});
   print OUT2 "$file2\_$file1\t$Di_nt\t$n_of_species\t$score{$Di_nt}\n";
   print "$file2\_$file1\t$Di_nt\t$n_of_species\t$score{$Di_nt}\n";

}
 
##************************************************************#
%score=();
%species=();

open OUT2, ">$file1.$file2.2.frequency";
#print OUT2 "guide piRNA from $file2\n";
#print OUT2 "target piRNA from $file1\n";
  %NTM=();
 `bowtie $file2 -r -a -v 1 -p 8 $file1.seq --suppress 1,4,6,7 | grep + > $file1.$file2.bowtie.out`;
  open IN, "$file1.$file2.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
   
   if ($_[3]=~/(\d):/) { next if ($1!=0);}
   $NTM{$_[2]}++;}

   open IN, "$file1.$file2.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
   if ($_[3] eq "") {$g_0_nt=substr($_[2],0,1); $t_9_nt=&revfa($g_0_nt); $score{$g_0_nt.$t_9_nt}+=$hash1{$file2}{$_[1]}*$hash2{$file1}{$_[2]}/$NTM{$_[2]};
                     $species{$g_0_nt.$t_9_nt}{$_[2]}=1;}
   elsif ($_[3]=~/(\d):(\w)>(\w)/) { next if ($1!=0);
     $t_9_nt=&revfa($2);$score{$3.$t_9_nt}+=$hash1{$file2}{$_[1]}*$hash2{$file1}{$_[2]}/$NTM{$_[2]};
     $species{$3.$t_9_nt}{$_[2]}=1;
     }}
   foreach $Di_nt(keys %score)
   {
      $n_of_species=scalar (keys %{$species{$Di_nt}});
      print OUT2 "$file1\_$file2\t$Di_nt\t$n_of_species\t$score{$Di_nt}\n";
      print "$file1\_$file2\t$Di_nt\t$n_of_species\t$score{$Di_nt}\n";
      }
 
     }
   

}
  #print "\n";
  #print OUT1 "\n";
  
}
 `rm *.ebwt`;
