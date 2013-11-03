#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
#p17-21 doesn't need to pair
#force 1mm at any position in p1-p16 0mm not included
# simplifized version with prefix 16nt
# input as norm.bed ( no header)
## separate each G1T10 combination (total 16)


open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while(<IN>) { if (/>(.+) type/) { $chr="chr$1";} else { chomp; $genome{$chr}=$_;}}


@dipairs=("AT","TA","GC","CG","AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");

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

 } }
}


# bowtie mapping and score calculating
for ($i=0; $i<$ARGV[2]; $i++)
{
   for ($j=0; $j<$ARGV[2]; $j++)   ###different
   {
         
      $file1=fileparse($ARGV[$i]); 
      @namefield=split(/\./,$file1);
      $name1=$namefield[2]."_".$namefield[1];
      $file1=$name1;
      $file2=fileparse($ARGV[$j]);
      @namefield=split(/\./,$file2);
      $name2=$namefield[2]."_".$namefield[1];
      $file2=$name2;
    
   open O, ">$file2.$file1.1mm.zscore.out";
   open OUT1,">$file2.$file1.1mm.pp_fraction";
      
   if ($total{$file1}<10 || $total{$file2}<10) { print "$file2-$file1\t-10"; print O "$file2-$file1\t-10"; print OUT1 "$file2-$file1\tNA";} 
   else {
   
         $X=0; $Z=0; %score=(); %species=();%pairs=();$count_N=0; %ppseq=();
          %X0=(); %Z0=();%s=(); %spe=();%pa=();%count_N0=(); 

        foreach ($n=1;$n<=20;$n++)
        {
        # file1 as ref
         `bowtie $file1.$n -r -a -v 1 -p 8 $file2.seq --suppress 1,4,6,7 | grep + > $file2.$file1.$n.bowtie.out`;
            foreach ($m=0;$m<=15;$m++)
            {
              %NTM=();
               open IN, "$file2.$file1.$n.bowtie.out";
               while(<IN>)
               {
                  chomp;
                  split(/\t/);  
                  if ($_[3]=~/(\d+):/)
                  {
                     if ($1==$m)
                     {$NTM{$_[2]}++;}
                  }  
               }
          
               open IN, "$file2.$file1.$n.bowtie.out";
               while(<IN>)
               {
                  chomp;
                  split(/\t/);
                  if ($_[3]=~/(\d+):(\w)>(\w)/)
                  {
                     if ($1==$m)
                     {
                     $score{$m}{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};
                     $species{$m}{$n}{$_[2]}=1;
                     $pairs{$m}{$n}++;
                     #$mm=$m+1;
                     #print OUT2 "$_[2]\tNA\t$mm\n" if ($n==10);
                     $ppseq{$m}{$_[2]}=1 if ($n==10);;
                     ###specific pairs (AU, GU...)
                     $t_9_nt=&revfa($2);
                     $s{$3.$t_9_nt}{$m}{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};
                     $spe{$3.$t_9_nt}{$m}{$n}{$_[2]}=1 ; ##number of species
                     $pa{$3.$t_9_nt}{$m}{$n}++; ## number of pairs
                     }
                  }
               }
            }
         }
        


         open OUT, ">$file2.$file1.1mm.pp";
         open OUT3, ">$file2.$file1.1mm.pp_reads_percentage";
         foreach ($m=0;$m<=15;$m++)
         {
            %count_N=();$X=0;$Z=0;$mm=$m+1;
            $std=0;$n_of_species=0;$n_of_s=0;$n_of_p=0;
            
            open OUT2, ">$file2.$file1.1mm.position_$mm.ppseq";
            foreach $seq (keys %{$ppseq{$m}})
            {
               print OUT2 "$seq\t$mm\n";
            }
            close(OUT2);
            
            foreach ($n=1;$n<=20;$n++)
            {         
            $score{$m}{$n}=0 if (!exists $score{$m}{$n});
            $n_of_species=scalar (keys %{$species{$m}{$n}});
            $n_of_pairs=$pairs{$m}{$n};
            print OUT "$n\t$n_of_species\t$n_of_pairs\t$score{$m}{$n}\t$mm\n";
            $count_N{$m}++ if ($score{$m}{$n}>0);
            }
         
        
           $X=$score{$m}{10}; delete $score{$m}{10};
           $n_of_s=scalar (keys %{$species{$m}{10}});
           $n_of_p=$pairs{$m}{10};
           $std=&std(values %{$score{$m}}); 
           if ($std>0 && $count_N{$m}>=5)
           
           { $Z=($X-&mean(values %{$score{$m}}))/$std;}
           
           else {$Z=-10;}
	   #$p=$m+1;
           print O "$file2-$file1\t$mm\t$n_of_s\t$n_of_p\t$Z\n";
           print "$file2-$file1\t$mm\t$n_of_s\t$n_of_p\t$Z\n";
         
           if ($Z!=-10) { print OUT1 "$file2-$file1\t$mm\t",$X/$total{$file1}/$total{$file2}*1000000,"\n";} else { print OUT1 "$file2-$file1\t$p\tNA\n";}
           $N1=`match.pl $file2.$file1.1mm.position_$mm.ppseq $file2.seq | sumcol+ 2`; chomp($N1);
           print OUT3 "mismatch at position","$mm\t", "$N1\t$total{$file2}\t",$N1/$total{$file2},"\n";
         }
         close (OUT3);
         close (OUT);
         close (OUT1);
         close (O);
         
         
      open OUT4, ">$file2.$file1.1U10A.1mm.zscore.out";
      open OUT5, ">$file2.$file1.1U10A.pp";
      foreach $p (@dipairs)
      {
         foreach ($m=0;$m<=15;$m++)
         {
            $X0=0;$Z0=0; $mm=$m+1;
            foreach ($n=1;$n<=20;$n++)
            {  
            ###specific pairs GU,AU...
            $s{$p}{$m}{$n}=0 if (!exists $s{$p}{$m}{$n});
            $nofs=scalar (keys %{$spe{$p}{$m}{$n}});
            $nofp=$pa{$p}{$m}{$n};
            print OUT5 "$n\t$nofs\t$nofp\t$s{$p}{$m}{$n}\t$mm\n";
            $count_N0{$p}{$m}++ if ($s{$p}{$m}{$n}>0);
            }
            $X0=$s{$p}{$m}{10}; delete $s{$p}{$m}{10};
            $std0=&std(values %{$s{$p}{$m}});
            if ($std0>0 && $count_N0{$p}>=5) { $Z0=($X0-&mean(values %{$s{$p}{$m}}))/$std0;} else {$Z0=-10;}
            $nof_s=scalar ( keys %{$spe{$p}{$m}{10}});
            $nof_p=$pa{$p}{$m}{10};
            print OUT4 "$file2-$file1\t$p\t$mm\t$nof_s\t$nof_p\t$Z0\n";  
         }   
      }
      close(OUT5);
      close(OUT4);     
   }

  }
}
`rm *.ebwt`;
