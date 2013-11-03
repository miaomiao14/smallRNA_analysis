#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
use POSIX;

##pp6
##input is xkxh.norm.bed instead of xkxh.transposon.mapper2

##output should include data for cross correlation analysis
## pp10 seq files:
## pp-16 seq files

## raw score

#09/11/2012

open IN, "/home/wangw1/pipeline/common/dm3.chrom.sizes";
while(<IN>)
{
   next if (/#/);
   chomp;split(/\t/);
   
   $chrsize{$_[0]}=$_[1];
   #for($i=0;$i<$_[1];$i++)
   #{
   #   ${$genome_pp10{$chr}}[$i]=0;
   #   ${$genome_pp16{$chr}}[$i]=0;
   #}
}


for ($i=0; $i<$ARGV[2]; $i++) {
   #print $ARGV[$i];

   for ($j=0; $j<$ARGV[2]; $j++) {
   $file1=fileparse($ARGV[$i]);  #$file1 is the target strand
   @namefield=split(/\./,$file1);
   $name1=$namefield[2]."_".$namefield[1];
   $file1=$name1;
    
   $file2=fileparse($ARGV[$j]); #file2 is the guide strand
   @namefield=split(/\./,$file2);
   $name2=$namefield[2]."_".$namefield[1];
   $file2=$name2;
    print "$file1-$file2";
   open PP1, ">$file1.$file2.product.pp";
     open PP2, ">$file1.$file2.min.pp";
      open PPSEQ, ">$file1.$file2.ppseq";
      open PPNEI, ">$file1.$file2.pp.neighbor.seq";
     $X=0; $Z=0; %score=();%s=();

     foreach ($n=-50;$n<=50;$n++) {
     %pp=(); %pos=(); %pp_seq=(); %pos_seq=();
     open IN, $ARGV[$i];
     while(<IN>) { chomp; split(/\t/);
      next if (/data/);
     #next if (length($_[0])>29 || length($_[0])<23); ##09/10/2012
     #$_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
     if ($_[3] eq '+') { $start=$_[1]+$n-1; $pp{"$_[0]:$start-"}+=$_[5]/$_[6];
                        #push @{$pp_seq{"$_[0]:$start-"}},[$_[0],$_[1],$_[3],$_[4]];
                        $pp_seq{"$_[0]:$start-"} = [$_[0],$_[1],$_[3],$_[4]]; #chr, start,strand,seq
                        }
     else { $start=$_[2]-$n+1;$pp{"$_[0]:$start+"}+=$_[5]/$_[6];
           #push @{$pp_seq{"$_[0]:$start+"}},[$_[0],$_[2],$_[3],$_[4]];
           $pp_seq{"$_[0]:$start+"}=[$_[0],$_[2],$_[3],$_[4]];
           }
     }
   
     open IN, $ARGV[$j];
     while(<IN>) { chomp; split(/\t/);
                  next if (/data/);
     #next if (length($_[0])>29 || length($_[0])<23); ##09/10/2012
     #$_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
     if ($_[3] eq '+') { $pos{"$_[0]:$_[1]+"}+=$_[5]/$_[6];
                        #push @{$pos_seq{"$_[0]:$_[1]+"}}, [$_[0],$_[1],$_[3],$_[4]];
                        $pos_seq{"$_[0]:$_[1]+"} = [$_[0],$_[1],$_[3],$_[4]];
                        }
     else { $pos{"$_[0]:$_[2]-"}+=$_[5]/$_[6];
           #push @{$pos_seq{"$_[0]:$_[2]-"}},[$_[0],$_[2],$_[3],$_[4]];
           $pos_seq{"$_[0]:$_[2]-"}=[$_[0],$_[2],$_[3],$_[4]]; #we just need to record the 5' coordinates, we don't care the 
           }
     }

     foreach (keys %pos)
     {
         if ($_ && exists $pp{$_})
         {
            $score{$n}+=$pos{$_}*$pp{$_} ;
            $s{$n}+=&min($pos{$_},$pp{$_});
            #$spaces=length($pos_seq{$_})-10;
            #printf PPSEQ "%$spaces";
            if ($n==10)
            {
               print PPSEQ "$pp_seq{$_}->[0]\t$pp_seq{$_}->[1]\t$pp_seq{$_}->[2]\t$pp_seq{$_}->[3]\t$pp{$_}\t$pos_seq{$_}->[0]\t$pos_seq{$_}->[1]\t$pos_seq{$_}->[2]\t$pos_seq{$_}->[3]\t$pos{$_}\n" ;
               # $pp_seq{$_}->[0]: the chromosome of guide strand
               # $pp_seq{$_}->[1]: the 5' end coordinate of guide strand
               # $pp_seq{$_}->[2]: the strand of guide strand
               # $pp_seq{$_}->[3]: the sequence of guide strand
               # $pp{$_} : the reads # of guide strand
               $pp10_mid=floor(($pp_seq{$_}->[1]+$pos_seq{$_}->[1])/2);
               #${$genome_pp10{$pp_seq{$_}->[0]}}[$pp10_mid]=1;
               ${$genome_pp10{$pp_seq{$_}->[0]}}{$pp10_mid}=1;
            }
            if ($n==-15)
            {
               print PPNEI "$pp_seq{$_}->[0]\t$pp_seq{$_}->[1]\t$pp_seq{$_}->[2]\t$pp_seq{$_}->[3]\t$pp{$_}\t$pos_seq{$_}->[0]\t$pos_seq{$_}->[1]\t$pos_seq{$_}->[2]\t$pos_seq{$_}->[3]\t$pos{$_}\n" ;
               # $pos_seq{$_}->[0]: the chromosome of guide strand
               # $pos_seq{$_}->[1]: the 5' end coordinate of guide strand
               # $pos_seq{$_}->[2]: the strand of guide strand
               # $pos_seq{$_}->[3]: the sequence of target strand
               # $pos{$_} : the reads # of target strand
               $pp16_mid=floor(($pp_seq{$_}->[1]+$pos_seq{$_}->[1])/2);
               #${$genome_pp16{$pp_seq{$_}->[0]}}[$pp16_mid]=1;
               ${$genome_pp16{$pp_seq{$_}->[0]}}{$pp16_mid}=1;
            }
         }
      }


     $score{$n}=0 if (!exists $score{$n});
     print PP1 "$n\t$score{$n}\n";
     $s{$n}=0 if (!exists $s{$n});
     print PP2 "$n\t$s{$n}\n";
     
     if ($n==10) { $X=$score{$n}; delete $score{$n};}
     }
   $std=&standard_deviation(values %score);
   if ($std>0) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
   print "\t$Z\t$X\n";
   close(PPNEI);
   close(PPSEQ);
   close (PP2);
   close (PP1);
   
     foreach $chr (keys %genome_pp10)
     {
      $q=0;
      
      open OUT_V, ">$file1.$file2.$chr.vectors_for_cross_correlation.out";
      
      print OUT_V "$chr\tpp10";
      
      foreach $p (sort {$a <=> $b} keys %{$genome_pp10{$chr}})
      {
         for ($i=$q;$i<$p;$i++)
         {
         print OUT_V "\t0";
         }
         print OUT_V "\t1";
         $q=$p+1;
      }
      
      for ($i=$q; $i< $chrsize{$chr}; $i++)
      {
         print OUT_V "\t0";
      }
      
      print OUT_V "\n";
      
      #$lastq=$q;
      $q=0;
      print OUT_V "$chr\tpp-16";
      
      $n_of_keys=scalar (%{$genome_pp16{$chr}});
      if($n_of_keys)
      {
         foreach $p (sort {$a <=> $b} keys %{$genome_pp16{$chr}})
         {
            for ($i=$q;$i<$p;$i++)
            {
            print OUT_V "\t0";
            }
            print OUT_V "\t1";
            $q=$p+1;
         }
         print OUT_V "\n";
      }

      for ($i=$q;$i<$chrsize{$chr};$i++)
      {
         print OUT_V "\t0";
      }   

      close(OUT_V);
     }
     
   
   }

}


sub mean {
my $count=0;
my(@numbers) =@_;
foreach (@_) { $count+=$_;}
return $count/(scalar @_);
}

sub standard_deviation {
my(@numbers) = @_;
#Prevent division by 0 error in case you get junk data
return undef unless(scalar(@numbers));

# Step 1, find the mean of the numbers
my $total1 = 0;
foreach my $num (@numbers) {
$total1 += $num;
}
my $mean1 = $total1 / (scalar @numbers);

# Step 2, find the mean of the squares of the differences
# between each number and the mean
my $total2 = 0;
foreach my $num (@numbers) {
$total2 += ($mean1-$num)**2;
}
my $mean2 = $total2 / (scalar @numbers);

# Step 3, standard deviation is the square root of the
# above mean
my $std_dev = sqrt($mean2);
return $std_dev;
}

