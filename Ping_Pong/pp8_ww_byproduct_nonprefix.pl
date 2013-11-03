#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
# rule is p1 and p17-21 doesn't need to pair but p2-10 need, the rest 1mm
# simplifized version with prefix 16nt
# input as norm.bed ( no header)

##check the byproducts of argonaute cleavage


open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while(<IN>) { if (/>(.+) type/) { $chr="chr$1";} else { chomp; $genome{$chr}=$_;}}

%hash1=();
%hash2=();

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
$hash2{$file}{$_[4]}+=$_[5]/$_[6];
for ($n=-50;$n<=50;$n++) {
   $len=length($_[4]);
if ($_[3] eq '+') { $start=$_[1]+$n; $str=substr($genome{$_[0]},$start,$len); $str=&revfa($str); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
else {  $start=$_[2]-$n; $str=substr($genome{$_[0]},$start,$len); $hash1{$file}{$n}{$str}+=$_[5]/$_[6];}
}
}
if ($total{$file}>10) {
open OUT, ">$file.seq";
foreach $seq (keys  %{$hash2{$file}}) { print OUT "$seq\t$hash2{$file}{$seq}\n" if($seq);}
close(OUT);
for ($n=-50;$n<=50;$n++) {
open OUT, ">$file.ref.$n.fa";
foreach (keys %{$hash1{$file}{$n}}) { print OUT ">$_\t$hash1{$file}{$n}{$_}\n$_\n" if($_);}
`bowtie-build $file.ref.$n.fa $file.$n`;
#`rm $file.ref.$n.fa`;
 }
close (OUT);
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

if ($total{$file1}<10 || $total{$file2}<10) { print "$file2-$file1\t-10"; print O "$file2-$file1\t-10"; print OUT1 "$file2-$file1\tNA\n";} 
else {

 $X=0; $Z=0; %score=();$count_N=0;
open OUT, ">$file2.$file1.pp";
open OUT2, ">$file2.$file1.ppseq";
foreach ($n=-50;$n<=50;$n++) {
# file1 as ref
 `bowtie $file1.$n -r -a -v 0 -p 8 $file2.seq --suppress 1,4,6,7 | grep + > $file2.$file1.$n.bowtie.out`;
  %NTM=();
  
  open IN, "$file2.$file1.$n.bowtie.out";
  while(<IN>) { chomp; split(/\t/);
  
     $score{$n}+=&min($hash1{$file1}{$n}{$_[1]},$hash2{$file2}{$_[2]});
     print OUT2 "$_[2]\n" if ($n==10);
  }  

     $score{$n}=0 if (!exists $score{$n});
     print OUT "$n\t$score{$n}\n";
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
