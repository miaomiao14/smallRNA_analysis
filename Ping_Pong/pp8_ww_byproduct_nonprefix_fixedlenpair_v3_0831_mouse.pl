#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
# rule is p1 and p17-21 doesn't need to pair but p2-10 need, the rest 1mm
# simplifized version with prefix 16nt
# input as norm.bed ( no header)

##check the byproducts of argonaute cleavage
##for Ago1
##no length restriction: 23-29nt


open IN, "/home/wangw1/pipeline/common/mmu_fasta/mm9_chromFa/mm9.fa.formated.fa";  ##the fasta file must be concatenated
while (<IN>)
{
    if(/>(.+)/)
    {
        $chr=$1;
    }
    else
    {
        chomp;
        $genome{$chr}=$_;
    }
}

%hash1=();
%hash2=();
%hash3=();

# make bowtie index
for ($i=0; $i<$ARGV[2]; $i++) {
my $file=fileparse($ARGV[$i]);
    @namefield=split(/\./,$file);
    $name=$namefield[2]."_".$namefield[1];
    push @argos, $name;
    $file=$name;
open IN, $ARGV[$i];
while(<IN>) { chomp; split(/\t/);  
#next if (length($_[4])>29 || length($_[4])<23);
next if (/data/);
$total{$file}+=$_[5]/$_[6];
$hash2{$file}{$_[4]}+=$_[5]/$_[6];
for ($n=-100;$n<=100;$n++) {
   $len=length($_[4]);
if ($_[3] eq '+') { $start=$_[1]+$n-($len+1); $str=substr($genome{$_[0]},$start,$len); $str=&revfa($str); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];$hash3{$file}{$n}{$str}=$_[4];}
else {  $start=$_[2]-$n; $str=substr($genome{$_[0]},$start,$len); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];$hash3{$file}{$n}{$str}=$_[4];}
}
}

if ($total{$file}>10) {
open OUT, ">$file.seq";
foreach (keys %{$hash2{$file}}) { print OUT "$_\t$hash2{$file}{$_}\n" ;}
close(OUT);}

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

if ($total{$file1}<10 || $total{$file2}<10) { print "$file2-$file1\t-10"; print O "$file2-$file1\t-10"; print OUT1 "$file2-$file1\tNA\n";} 
else {

 $X=0; $Z=0; %s_mul=();$count_N=0; %s_min=();%s_max=();%s_mean=();%score=();
open OUT5, ">$file2.$file1.min.pp";
open OUT6, ">$file2.$file1.max.pp";
open OUT7, ">$file2.$file1.mean.pp";
open OUT8, ">$file2.$file1.mul.pp";

open OUT2, ">$file2.$file1.ppseq";
open OUT21, ">$file2.$file1.16merbyproduct.ppseq";
foreach ($n=-100;$n<=100;$n++) {
  %NTM=();
  foreach $seq (keys %{$hash1{$file2}{$n}}) ##file2 guide
  {
   
      if($hash2{$file1}{$seq}) ## file1 target
      {
      $s_min{$n}+=&min($hash1{$file2}{$n}{$seq},$hash2{$file1}{$seq});
      $s_max{$n}+=&max($hash1{$file2}{$n}{$seq},$hash2{$file1}{$seq});
      $s_mean{$n}+=&mean($hash1{$file2}{$n}{$seq},$hash2{$file1}{$seq});
      $s_mul{$n}+=$hash1{$file2}{$n}{$seq}*$hash2{$file1}{$seq};
      $score{$n}+=$hash1{$file2}{$n}{$seq}*$hash2{$file1}{$seq};
      print OUT2 "$seq\t$hash3{$file2}{$n}{$str}\n" if ($n==10);
      print OUT21 "$seq\t$hash3{$file2}{$n}{$str}\n"if ($n==-16);
      }
  }  

     $score{$n}=0 if (!exists $score{$n});
     print OUT5 "$n\t$s_min{$n}\n";
     print OUT6 "$n\t$s_max{$n}\n";
     print OUT7 "$n\t$s_mean{$n}\n";
     print OUT8 "$n\t$s_mul{$n}\n";
     $count_N++ if ($score{$n}>0);
}
   $X=$score{10}; delete $score{10};
   $std=&std(values %score); 
   if ($std>0 && $count_N>=5) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
   print O "$file2-$file1\t\t$Z\n";
   print "$file2-$file1\t\t$Z\n";
 
   if ($Z!=-10) { print OUT1 "\t",$X/$total{$file1}/$total{$file2}*1000000,"\n";} else { print OUT1 "\tNA\n";}
   $N1=`match.pl $file2.$file1.ppseq $file2.seq | sumcol+ 2`; chomp($N1);
   open OUT3, ">$file2.$file1.pp_reads_percentage";
   print OUT3 "$N1\t$total{$file2}\t",$N1/$total{$file2},"\n";
}

}

}
`rm *.ebwt`;
