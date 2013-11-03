#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
# rule is p1 doesn't need to pair
# 5'end overlap 1,2,..,9,11,...20 as backgroud noise
# input as norm.bed ( no header)


open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while(<IN>) { if (/>(.+) type/) { $chr="chr$1";} else { chomp; $genome{$chr}=$_;}}
@pairs=("AT","TA","GC","CG","AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");

# make bowtie index
for ($i=0; $i<@ARGV; $i++) {
my $file=fileparse($ARGV[$i]);
    @namefield=split(/\./,$file);
    $name=$namefield[2]."_".$namefield[1];
    push @argos, $name;
    $file=$name;
#%hash1=();
#%hash2=();
#%total=();
    
open IN, $ARGV[$i];
while(<IN>) { chomp; split(/\t/);  
next if (length($_[4])>29 || length($_[4])<23);
next if (/data/);
$total{$file}+=$_[5]/$_[6];
#$hash2{$file}{substr($_[4],0,16)}+=$_[5]/$_[6];
for ($n=1;$n<=20;$n++) {
if ($_[3] eq '+')
    {
        #$start=$_[1]+$n-1;
        $start=$_[1]-1;
        $str=substr($genome{$_[0]},$start,$n); $str=&revfa($str);
        $hash1{$file}{$n}{$str}+=$_[5]/$_[6];
        $hash2{$file}{$n}{substr($_[4],0,$n)}+=$_[5]/$_[6];
    }
else
    {
        $start=$_[2]-$n;
        $str=substr($genome{$_[0]},$start,$n);
        $hash1{$file}{$n}{$str}+=$_[5]/$_[6];
        $hash2{$file}{$n}{substr($_[4],0,$n)}+=$_[5]/$_[6];
        }
}
}
if ($total{$file}>10) {

for ($n=1;$n<=20;$n++) {

open OUT, ">$file.$n.seq";
foreach (keys %{$hash2{$file}{$n}}) { print OUT "$_\t$hash2{$file}{$n}{$_}\n" if(length($_)==$n);}    
    
open OUT, ">$file.ref.$n.fa";
foreach (keys %{$hash1{$file}{$n}}) { print OUT ">$_\t$hash1{$file}{$n}{$_}\n$_\n" if(length($_)==$n);}
`bowtie-build $file.ref.$n.fa $file.$n`;
#`rm $file.ref.$n.fa`;
 } }
}

#####################################################################################
open OUT1,">pp_fraction";
print OUT1 "pp";
print OUT1 map {"\t$_"} @argos;
print OUT1 "\n";
open O, ">zscore.out";
open OUA, ">UA_VA.zscore.out";
# bowtie mapping and score calculating
for ($i=0; $i<@ARGV; $i++) {
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
   #print OUA "$file1-$file2";
   print "$file1-$file2";
   print OUT1 "$file1-$file2";
   
if ($total{$file1}<10 || $total{$file2}<10) { print "\t-10"; print O "\t-10"; print OUA "$file1-$file2\t-10\n"; print OUT1 "\tNA";} 
else {

 $X=0; $Z=0; %score=();$count_N=0;
 %X0=(); %Z0=(); %s=(); %suv=(); %count_N0=();
open OUT, ">$file1.$file2.pp";
open OUT0, ">$file1.$file2.VA.pp";
open OUT2, ">$file1.$file2.ppseq";
foreach ($n=1;$n<=20;$n++)
{
# file1 as ref
 `bowtie $file1.$n -r -a -v 1 -p 8 $file2.$n.seq --suppress 1,4,6,7 | grep + > $file1.$file2.$n.bowtie.out`;
  %NTM=();
  open IN, "$file1.$file2.$n.bowtie.out";
  while(<IN>)
  {
    chomp; split(/\t/);
    
    if ($_[3]=~/(\d+):/) { next if ($1!=($n-1));}  # allow 1mm at the 10th position of target strand
    $NTM{$_[2]}++;
    #$count_N++;
    }
  open IN, "$file1.$file2.$n.bowtie.out";
  while(<IN>)
  {
    chomp; split(/\t/);
    if ($_[3] eq "")
    {
    $m=$n-1;
    $t_9_nt=substr($_[2],$m,1); $g_0_nt=&revfa($t_9_nt);
    $s{$g_0_nt.$t_9_nt}{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$n}{$_[2]}/$NTM{$_[2]};
    $score{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$n}{$_[2]}/$NTM{$_[2]};
    print OUT2 "NA\t$_[2]\n" if ($n==10);
    }
    elsif ($_[3]=~/(\d):(\w)>(\w)/)
    {
    next if ($1!=($n-1));  # allow 1mm at the 10th position of target strand
    $g_0_nt=&revfa($2);
    $s{$g_0_nt.$3}{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$n}{$_[2]}/$NTM{$_[2]};
    #$suv{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$n}{$_[2]}/$NTM{$_[2]};
    $score{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$n}{$_[2]}/$NTM{$_[2]};
    print OUT2 "NA\t$_[2]\n" if ($n==10);
    }
  }
 
# `rm $file1.$file2.$n.bowtie.out`;  
# file2 as ref
  %NTM=();
 `bowtie $file2.$n -r -a -v 1 -p 8 $file1.$n.seq --suppress 1,4,6,7 | grep + > $file2.$file1.$n.bowtie.out`;
 open IN, "$file2.$file1.$n.bowtie.out";
  while(<IN>)
  {
    chomp; split(/\t/);  
    if ($_[3]=~/(\d+):/) { next if ($1!=($n-1));}  # allow 1mm at the 10th position of target strand
    $NTM{$_[2]}++;
    # $count_N++;
  }
   open IN, "$file2.$file1.$n.bowtie.out";
  while(<IN>)
  {
    chomp; split(/\t/);
  
    if ($_[3] eq "")
    {
    $m=$n-1;
    $t_9_nt=substr($_[2],$m,1); $g_0_nt=&revfa($t_9_nt);
    $s{$g_0_nt.$t_9_nt}{$n}+=$hash1{$file2}{$n}{$_[1]}*$hash2{$file1}{$n}{$_[2]}/$NTM{$_[2]};
    $score{$n}+=$hash1{$file2}{$n}{$_[1]}*$hash2{$file1}{$n}{$_[2]}/$NTM{$_[2]};
    print OUT2 "NA\t$_[2]\n" if ($n==10);
    }
    elsif ($_[3]=~/(\d):(\w)>(\w)/)
    {
    next if ($1!=($n-1));  # allow 1mm at the 10th position of target strand
    $g_0_nt=&revfa($2);
    $s{$g_0_nt.$3}{$n}+=$hash1{$file2}{$n}{$_[1]}*$hash2{$file1}{$n}{$_[2]}/$NTM{$_[2]};
    #$suv{$n}+=$hash1{$file2}{$n}{$_[1]}*$hash2{$file1}{$n}{$_[2]}/$NTM{$_[2]};
    $score{$n}+=$hash1{$file2}{$n}{$_[1]}*$hash2{$file1}{$n}{$_[2]}/$NTM{$_[2]};
    print OUT2 "NA\t$_[2]\n" if ($n==10);
    }
# `rm $file2.$file1.$n.bowtie.out`;
  }

     $score{$n}=0 if (!exists $score{$n});
     
     foreach $p (@pairs)
     {
     $s{$p}{$n}=0 if (!exists $s{$p}{$n});
     print OUT0 "$n\t$p\t$s{$p}{$n}\n";
     $count_N0{$p}++ if ($s{$p}{$n}>0);
     }
     
     print OUT "$n\t$score{$n}\n";
     $count_N++ if ($score{$n}>0);

}  
   $X=$score{10}; delete $score{10};
   $std=&std(values %score); 
   if ($std>0 && $count_N>=5) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
   print O "\t$Z";
   print "\t$Z";
   
    foreach $p (@pairs)
    {
    #$std=0;
    $X0{$p}=$s{$p}{10}; delete $s{$p}{10};
    $std0=&std(values %{$s{$p}});
    if ($std0>0 && $count_N0{$p}>=5) { $Z0{$p}=($X0{$p}-&mean(values %{$s{$p}}))/$std0;} else {$Z0{$p}=-10;}
    print OUA "$file1-$file2\t$p\t$Z0{$p}\n";
    
    }
 
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
