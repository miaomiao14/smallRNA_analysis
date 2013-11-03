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
#$hash2{$file}{substr($_[4],0,16)}+=$_[5]/$_[6];
$hash2{$file}{$_[7]}{substr($_[4],0,16)}+=$_[5]/$_[6];
for ($n=1;$n<=20;$n++) {
if ($_[3] eq '+') { $start=$_[1]+$n-17; $str=substr($genome{$_[0]},$start,16); $str=&revfa($str); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
else {  $start=$_[2]-$n; $str=substr($genome{$_[0]},$start,16); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
}
}
if ($total{$file}>10) {
open OUT, ">$file.seq";
foreach $sense (keys %{$hash2{$file}})
{
   foreach $seq (keys %{$hash2{$file}{$sense}})
   {
   print OUT "$seq\t$hash2{$file}{$sense}{$seq}\t$sense\n" if (length($seq)==16);
   }
}
close(OUT);

open OUT, ">$file.seq.fa";
foreach $sense (keys %{$hash2{$file}})
{
   foreach $seq (keys %{$hash2{$file}{$sense}})
   {
   print OUT ">$sense\t$hash2{$file}{$sense}{$seq}\n$seq\n" if (length($seq)==16);
   }
}
   for ($n=1;$n<=20;$n++)
   {
      open OUT, ">$file.ref.$n.fa";
      foreach (keys %{$hash1{$file}{$n}})
      {
         print OUT ">$_\t$hash1{$file}{$n}{$_}\n$_\n" if (length($_)==16);$k++;
      }
      `bowtie-build $file.ref.$n.fa $file.$n`;
      #`rm $file.ref.$n.fa`;
      close (OUT);
   }

}
}


# bowtie mapping and score calculating
for ($i=0; $i<$ARGV[2]; $i++) {
   for ($j=0; $j<$ARGV[2]; $j++) {
      
   $file1=fileparse($ARGV[$i]); 
   @namefield=split(/\./,$file1);
   $name1=$namefield[2]."_".$namefield[1];
   $file1=$name1;
   $file2=fileparse($ARGV[$j]);
   @namefield=split(/\./,$file2);
   $name2=$namefield[2]."_".$namefield[1];
   $file2=$name2;
   

   
   open OUT1,">pp_fraction";
   open O, ">zscore.out";

if ($total{$file1}<10 || $total{$file2}<10) {  print O "$file2-$file1\t-10"; print "$file2-$file1\tNA\tNA\tNA\tNA\tNA\tNA\n"; print OUT1 "$file2-$file1\tNA\n";} 
else {

 $X=0; $Z=0; %score=();$count_N=0;
open OUT, ">$file2.$file1.pp";
open OUT2, ">$file2.$file1.ppseq";
foreach ($n=1;$n<=20;$n++) {
# file1 as ref
 `bowtie $file1.$n -a -v 1 -p 8 -f $file2.seq.fa --suppress 4,6,7 | grep + > $file2.$file1.$n.bowtie.out`;
  %NTM=();
  open IN, "$file2.$file1.$n.bowtie.out";
  while(<IN>) { chomp; split(/\t/);  
  if ($_[4]=~/(\d+):/) { next if ($1<=9 && $1>=1);}  # no seed mismatches 0 based
  $NTM{$_[0]}{$_[3]}++;
  #$count_N++;
  }
  open IN, "$file2.$file1.$n.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
  if ($_[4]=~/(\d+):/) { next if ($1<=9 && $1>=1);}  # no seed mismatches 0 based
     $score{$n}+=$hash1{$file1}{$n}{$_[2]}*$hash2{$file2}{$_[0]}{$_[3]}/$NTM{$_[0]}{$_[3]};
     print OUT2 "$_[3]\n" if ($n==10);
     $ppseq{$_[0]}{$_[3]}=1 if ($n==10);
  }  

     $score{$n}=0 if (!exists $score{$n});
     print OUT "$n\t$score{$n}\n";
     $count_N++ if ($score{$n}>0);
}
   $X=$score{10}; delete $score{10};
   $std=&std(values %score); 
   if ($std>0 && $count_N>=5) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
   print O "$file2-$file1\t$Z\n";
   #print "\t$Z\n";
 
   if ($Z!=-10) { print OUT1 "\t",$X/$total{$file1}/$total{$file2}*1000000,"\n";} else { print OUT1 "\tNA\n";}
   
   
   
   foreach $sense (keys %ppseq)
   {
      foreach $seq (keys %{$ppseq{$sense}})
      {
         $ppfraction{$sense}+=$hash2{$file2}{$sense}{$seq}; ##pp sense and antisense
      }
   }
   foreach $sense (keys %{$hash2{$file2}})
   {
      foreach $seq (keys %{$hash2{$file2}{$sense}})
      {
         if(! exists $ppseq{$sense}{$seq}) ##nonpp sense and antisense
         {
            $nonppfraction{$sense}+=$hash2{$file2}{$sense}{$seq};
         }
         $totalfraction{$sense}+=$hash2{$file2}{$sense}{$seq};
      }
   }
   if ($ppfraction{sense}==0 || $ppfraction{antisense}==0)
   {
      $pp_fra="NA";
   }
   else
   {
   $pp_fra=$ppfraction{sense}/($ppfraction{sense}+$ppfraction{antisense});
   }
   $nonpp_fra=$nonppfraction{sense}/($nonppfraction{sense}+$nonppfraction{antisense});
   $total_fra=$totalfraction{sense}/($totalfraction{sense}+$totalfraction{antisense});
   
   $N1=$ppfraction{sense}+$ppfraction{antisense};
   
   open OUT3, ">$file2.$file1.pp_reads_percentage";
   print OUT3 "$N1\t$total{$file2}\t",$N1/$total{$file2},"\n";
   
   print "$file2-$file1\t$pp_fra\t$nonpp_fra\t$total_fra\t","$N1\t$total{$file2}\t",$N1/$total{$file2},"\n";

}

}

}
`rm *.ebwt`;
