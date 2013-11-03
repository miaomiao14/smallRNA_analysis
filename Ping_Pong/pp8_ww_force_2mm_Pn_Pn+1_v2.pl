#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
#p17-21 doesn't need to pair
#force 2 consecutive mms at any position in p1-p16 0mm not included
# simplifized version with prefix 16nt
# input as norm.bed ( no header)
# separate Ago3:Aub and Aub:Ago3 PP
# record the number of PP pairs
#08/02/2012


open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while(<IN>) { if (/>(.+) type/) { $chr="chr$1";} else { chomp; $genome{$chr}=$_;}}
#$m=$ARGV[3];
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


# bowtie mapping and score calculating
for ($i=0; $i<$ARGV[2]; $i++)
{
   for ($j=0; $j<=$i; $j++)
   {
         
      $file1=fileparse($ARGV[$i]); 
      @namefield=split(/\./,$file1);
      $name1=$namefield[2]."_".$namefield[1];
      $file1=$name1;
      $file2=fileparse($ARGV[$j]);
      @namefield=split(/\./,$file2);
      $name2=$namefield[2]."_".$namefield[1];
      $file2=$name2;
      
      #print O "$file1-$file2";
      #print "$file1-$file2";
      #print OUT1 "$file1-$file2";
      
   if ($total{$file1}<10 || $total{$file2}<10) { print "$file1-$file2\t-10"; print O "$file1-$file2\t-10"; print OUT1 "$file1-$file2\tNA";} 
   else {
   
         $X=0; $Z=0; %score=(); %species=();
         open OUT, ">$file1.$file2.1.2cmm.pp";
         open OUT2, ">$file1.$file2.1.2cmm.ppseq";
         open OUT1,">$file1.$file2.1.2cmm.pp_fraction";
         #print OUT1 "pp";
         open O, ">$file1.$file2.1.2cmm.zscore.out";
         
        foreach ($n=1;$n<=20;$n++)
        {
        # file1 as ref
         `bowtie $file1.$n -r -a -v 2 -p 8 $file2.seq --suppress 1,4,6,7 | grep + > $file1.$file2.$n.bowtie.out`;
            foreach ($m=0;$m<=14;$m++)
            {
              %NTM=();
               open IN, "$file1.$file2.$n.bowtie.out";
               while(<IN>)
               { chomp; split(/\t/);  
                  if ($_[3]=~/(\d+):.*,(\d+)/)
                  {
                     if ($1==$m && $2==$m+1)
                     {$NTM{$_[2]}++;}
                  }  
               }
          
               open IN, "$file1.$file2.$n.bowtie.out";
               while(<IN>)
               {
               chomp; split(/\t/);
                  if ($_[3]=~/(\d+):.*,(\d+)/)
                  {
                     if ($1==$m && $2==$m+1)
                     {
                     $score{$m}{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};
                     $species{$m}{$n}{$_[2]}=1 ;
                     $p=$m+1;
                     print OUT2 "NA\t$_[2]\t$p\n" if ($n==10);
                     }
                  }
               }
            } #m
        } #n
        foreach ($m=0;$m<=14;$m++)
         {
            %count_N=();$X=0;$Z=0;$p=$m+1;$q=$p+1;
            $std=0;
            foreach ($n=1;$n<=20;$n++)
            {         
            $score{$m}{$n}=0 if (!exists $score{$m}{$n});
            $n_of_species=scalar (keys %{$species{$m}{$n}});
            print OUT "$n\t$n_of_species\t$score{$m}{$n}\t$p\t$q\n";
            $count_N{$m}++ if ($score{$m}{$n}>0);
            }
         
        
           $X=$score{$m}{10}; delete $score{$m}{10};
           $n_of_s=scalar (keys %{$species{$m}{10}});
           $std=&std(values %{$score{$m}}); 
           if ($std>0 && $count_N{$m}>=5) { $Z=($X-&mean(values %{$score{$m}}))/$std;} else {$Z=-10;}
	   #$p=$m+1;
           print O "$file2-$file1\t$p\t$q\t$n_of_s\t$Z\n";
           print "$file2-$file1\t$p\t$q\t$n_of_s\t$Z\n";
         
           if ($Z!=-10) { print OUT1 "$file2-$file1\t$p\t",$X/$total{$file1}/$total{$file2}*1000000,"\n";} else { print OUT1 "$file2-$file1\t$p\tNA\n";}
           #$N1=`match.pl $file1.$file2.1.2cmm.ppseq $file1.seq | sumcol+ 2`; chomp($N1);
           `cut -f2 $file1.$file2.1.2cmm.ppseq > $file1.$file2.1.2cmm.$m.ppseq.temp`;
           $N2=`match.pl $file1.$file2.1.2cmm.$m.ppseq.temp $file2.seq | sumcol+ 2`; chomp($N2);
           open OUT3, ">$file1.$file2.1.2cmm.$m.pp_reads_percentage";
           print OUT3 "mismatch","$p\t$q\t", "$N2\t$total{$file2}\t",$N2/$total{$file2},"\n";

         }


      $X=0; $Z=0; %score=(); %species=();
      open OUT, ">$file2.$file1.2.2cmm.pp";
      open OUT2, ">$file2.$file1.2.2cmm.ppseq";
      open OUT1,">$file2.$file1.2.2cmm.pp_fraction";
         #print OUT1 "pp";
      open O, ">$file2.$file1.2.2cmm.zscore.out";
      foreach ($n=1;$n<=20;$n++)
      {
        # `rm $file1.$file2.$n.bowtie.out`;  
        # file2 as ref
      
         `bowtie $file2.$n -r -a -v 2 -p 8 $file1.seq --suppress 1,4,6,7 | grep + > $file2.$file1.$n.bowtie.out`;
         foreach ($m=0;$m<=14;$m++)
         {
            %NTM=();
            open IN, "$file2.$file1.$n.bowtie.out";
            while(<IN>)
            {
               chomp; split(/\t/);  
               if ($_[3]=~/(\d+):.*,(\d+)/)
               {
                  if ($1==$m && $2==$m+1)
                  {
                  $NTM{$_[2]}++;
                  }
               }
            }
            open IN, "$file2.$file1.$n.bowtie.out";
            while(<IN>)
            { chomp; split(/\t/);
               if ($_[3]=~/(\d+):.*,(\d+)/)
               {
                  if ($1==$m && $2==$m+1)
                  {
                  $score{$m}{$n}+=$hash1{$file2}{$n}{$_[1]}*$hash2{$file1}{$_[2]}/$NTM{$_[2]};
                  $species{$m}{$n}{$_[2]}=1 ;
                  $p=$m+1;
                  print OUT2 "$_[2]\tNA\t$p\n" if ($n==10);
                  }
               }
            }
        # `rm $file2.$file1.$n.bowtie.out`;

         } #m
      } #n
        foreach ($m=0;$m<=14;$m++)
         {
            %count_N=();$X=0;$Z=0;$p=$m+1;$q=$p+1;
            $std=0;
            foreach ($n=1;$n<=20;$n++)
            {         
            $score{$m}{$n}=0 if (!exists $score{$m}{$n});
            $n_of_species=scalar (keys %{$species{$m}{$n}});
            print OUT "$n\t$n_of_species\t$score{$m}{$n}\t$p\t$q\n";
            $count_N{$m}++ if ($score{$m}{$n}>0);
            }
         
        
           $X=$score{$m}{10}; delete $score{$m}{10};
           $n_of_s=scalar (keys %{$species{$m}{10}});
           $std=&std(values %{$score{$m}}); 
           if ($std>0 && $count_N{$m}>=5) { $Z=($X-&mean(values %{$score{$m}}))/$std;} else {$Z=-10;}
	   #$p=$m+1;
           print O "$file1-$file2\t$p\t$q\t$n_of_s\t$Z\n";
           print "$file1-$file2\t$p\t$q\t$n_of_s\t$Z\n";
         
           if ($Z!=-10) { print OUT1 "$file1-$file2\t$p\t",$X/$total{$file1}/$total{$file2}*1000000,"\n";} else { print OUT1 "$file1-$file2\t$p\tNA\n";}
           $N1=`match.pl $file2.$file1.2.2cmm.ppseq $file1.seq | sumcol+ 2`; chomp($N1);
           #`cut -f2 $file2.$file1.2.2cmm.ppseq > $file2.$file1.2.2cmm.$m.ppseq.temp`;
           #$N2=`match.pl $file2.$file1.2.2cmm.$m.ppseq.temp $file2.seq | sumcol+ 2`; chomp($N2);
           open OUT3, ">$file2.$file1.2.2cmm.$m.pp_reads_percentage";
           print OUT3 "mismatch","$p\t$q\t", "$N1\t$total{$file1}\t",$N1/$total{$file1},"\n";

         }
        }

        }
   }
`rm *.ebwt`;
