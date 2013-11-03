#!/usr/bin/perl
use File::Basename;
BEGIN { unshift @INC,"/home/xuj1/bin";}
require "Statstics.pm";
require "Jia.pm";



##report PP pair sequences
##WEI WANG
##09/2012

## when report mean and std, need to consider sequencing depth

##09/25/2012
# For N1, if it can find M10, (their 5' end distance is 10 nt apart), then for d=1,2..9,11..20,
# if d=1, for N1, min(max(N1M1-N1M10,0),max(0,N1M1-M1N10))
# if d=2, for N1, min(max(N1M2-N1M10,0), max(0,N1M2-M2N11)) if exists M2N11
# if d=3, for N1, min(max(N1M3-N1M10,0),max(0,N1M3-M3N12)) if exists M3N12

##v3##
##09/27/2012
## Phil is out today
## need to take the noise-only species into consideration

##p9##
##09/29/2012
##BO's idea, check if the refunded part is artefact or not, refund p9, p11...
##How can we say p10 is signal while others are noises?


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
   
   open PPSEQ, ">$file1_$file2.refund.ppseq";
   open PP, ">$file1_$file2.refund.pp";

   open PPZ, ">$file1_$file2.refund.zscore.out";

   
   print PPZ "$file1-$file2\t";
   
   $X=0; $Z=0; %score=(); %s=();


   %pp=(); %pos=(); %pp_seq=(); %pos_seq=();
   %pp_10=();
      open IN, $ARGV[$i];  ##suppose Ago3IP
      while(<IN>)
      {
         chomp;
         s/\s+/\t/g;
         split(/\t/);
         next if (length($_[0])>29 || length($_[0])<23);
         $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
         if ($4 eq '+')
         {
            foreach ($n=1;$n<=32;$n++)
            {
               $start=$2+$n-1;
               $pp{"$1:$2+"}{$n}{"$1:$start-"}+=$_[1]/$_[6]; ##record each species, and their corresponding pairs when 5' end distance is 1,2,...20
                                                            ## the reads number is from N1
               
               $pp_10{"$1:$2+"}="$1:$start-" if($n==11);  #record pp9 for decision, key is the N1, value is the M1
               #$pp_seq{"$1:$start-"}=$_[0] if($n==10); ## to make sure the key are unique, won't be overwrite by some non-10 nt overlap sequence pair
            } #shift n
         }
         else
         {
            foreach ($n=1;$n<=32;$n++)
            {
               $start=$3-$n+1;
               $pp{"$1:$3-"}{$n}{"$1:$start+"}+=$_[1]/$_[6]; ##record each species, and their corresponding pairs when 5' end distance is 1,2,...20
                                                            ## the reads number is from N1
               $pp_10{"$1:$3-"}="$1:$start+" if($n==11); #record pp9 for decision
               #$pp_seq{"$1:$start+"}=$_[0] if($n==10); ## to make sure the key are unique, won't be overwrite by some non-10 nt overlap sequence pair
            } #shift n
         }
         #print "hello\n";
      }
   
      open IN, $ARGV[$j]; ##suppose AubIP
      while(<IN>)
      {
         chomp; s/\s+/\t/g;split(/\t/);
         next if (length($_[0])>29 || length($_[0])<23);
         $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
         if ($4 eq '+')                   
         {
            $pos{"$1:$2+"}+=$_[1]/$_[6]; # original species
            #$pos_seq{"$1:$2+"}= $_[0];
         }
         else
         {
            $pos{"$1:$3-"}+=$_[1]/$_[6]; # original species
            #$pos_seq{"$1:$3-"}=$_[0];
         }
      }
####################################################################################################################################################################################
  
   foreach $species (keys %pp) ##For N1 (the key), if it can find M10, (their 5' end distance is 10 nt apart), then N1M2, N1M3 should be refunded!
   ## the key of %pp_10 is the same as the key of the first hash layer of %pp
   {
         
         ##if we only consider those species which have p10 pair??
         ## in Jia's pp6, no matter the p10 exists or not, the background noise will be counted.
         ## Here, I only take the noise around pp10 into account
         
         ###adjusted in v3
            foreach $n (keys %{$pp{$species}})
            {  

                  foreach $species_pair (keys %{$pp{$species}{$n}}) ##M1,M2,...M20
                  {
                     $temp1=0;
                     $temp2=0;
                     if($pos{$species_pair}) #if find M1, 
                     {
                        if($n!=11 && $pos{$pp_10{$species}}) #if n=10, no need to refund
                        {
                           #if($pos{$pp_10{$species}}) #find a 10 nt pair M10 (value is M10)
                           #{
                              $temp1=&max(0,$pp{$species}{$n}{$species_pair}*($pos{$species_pair}-$pos{$pp_10{$species}}) ); #N1M1-N1M10 for n=1
                           #}
                        }
                        else
                        {
                        #############################################################################
                        #########################For test, for ping-pong species, times 10############
                        ##############################################################################
                        $temp1=$pp{$species}{$n}{$species_pair}*$pos{$species_pair}; #$temp1=N1M1  ##if $pos{$species_pair} does not exist, it's 0.
                        }
                        ##if($pos{$pp_10{$species_pair}} && $n!=10) ##find a 10nt pair for M1,M2,...M20
                        ##$pos{$pp_10{$species_pair}} this is not looking for pp10 pairs from the other IP datasets, but from the same one
                        
                        ##Hopefully, the following is the right one
                        $temp2=$pp{$species}{$n}{$species_pair}*$pos{$species_pair};
                        if ($pp{$pp_10{$species_pair}}{11}{$species_pair}) #for M1 if N10 exists
                        
                        {
                           if ($n!=11) #if n=10, no need to refund
                           {
                           $temp2=&max(0,$pos{$species_pair}*($pp{$species}{$n}{$species_pair}-$pp{$pp_10{$species_pair}}{11}{$species_pair}) ); #M1N1-M1N10 for n=1
                           }
                           else
                           {
                           ##############################################################################
                           #########################For test, for ping-pong species, times 10############
                           ##############################################################################
                           $temp2=$pp{$species}{$n}{$species_pair}*$pos{$species_pair}; #$temp2=N1M1
                           }
                        }
                     $score{$n}+=&max($temp1,$temp2);
                     #print PPSEQ "$pp_seq{$species_pair}\t$pp{$species}{$n}{$species_pair}\t$pos_seq{$species_pair}\t$pos{$species_pair}\n" if ($n==10);
                     }    
                     
                  }
                  

               
            }
            
         #}#  
   }
####################################################################################################################################################################################     
  
   
   foreach ($n=1;$n<=32;$n++)
   {
      
      $score{$n}=0 if (!exists $score{$n});
      print PP "$n\t$score{$n}\n";
      #if ($n==10) { $X=$score{$n}; delete $score{$n};}
   }
   
      $X=$score{10};
      $p9=$score{9};
      $p11=$score{11};
      
      %score_20=();
      for($j=1;$j<=20;$j++)
      {
         $score_20{$j}=$score{$j};
      }
      
      
      $total_score=0;
      $total_score=&total(values %score);
      if ($total_score!=0)
      {
         $percentage=$p11/$total_score;
      }
      else
      {
         $percentage=0;
      }
      delete $score{11};
      $std=&standard_deviation(values %score);
      $m=&mean(values %score);
      if($p11>1)  ##if the signal itself is very low
      {
         if ($std>0) ##in terms of refund, when signals are very pure, std=0
         {
            $Z=($p11-&mean(values %score))/$std;
         }
         else
         {
            $Z=10000;
         }
      }
      else
      {
         $Z=-10;
      }
      
      ##how to normalzie p10,p11,p9, mean and std
      ## because they are the sum of multiplification
      
      print PPZ "$Z\t$percentage\t$m\t$std\t$X\t$p9\t$p11\t";
      print "$Z\t$percentage\t$m\t$std\t$X\t$p9\t$p11\t";
      
      $total_score=0;
      $total_score=&total(values %score_20);
      if ($total_score!=0)
      {
         $percentage=$p11/$total_score;
      }
      else
      {
         $percentage=0;
      }
      delete $score_20{11};
      $std=&standard_deviation(values %score_20);
      $m=&mean(values %score_20);
      
      if($p11>1)  ##if the signal itself is very low
      {
         if ($std>0)  ##in terms of refund, when signals are very pure, std=0
         {
            $Z=($p11-&mean(values %score_20))/$std;
         }
         else
         {
            $Z=10000;
         }
      }
      else
      {
         $Z=-10;
      }
      print PPZ "$Z\t$percentage\t$m\t$std\n";
      print "$Z\t$percentage\t$m\t$std\n";
      
   
   
   #foreach $species (keys %s)
   #{
   #   foreach ($n=1;$n<=20;$n++)
   #   {
   #      $s{$species}{$n}=0 if (!exists $s{$species}{$n});
   #      print PPlocal "$species\t$n\t$s{$species}{$n}\n";
   #      if ($n==10) { $X=$s{$species}{$n}; delete $s{$species}{$n};}
   #   }
   #   $std=&standard_deviation(values %{$s{$species}});
   #   
   #   if ($std>0)
   #   {
   #      $Z=($X-&mean(values %{$s{$species}}))/$std;
   #   }
   #   else
   #   {
   #      $Z=-10;
   #   }
   #   print PPZlocal "$species\t$Z\t$X\t$std\n";   
   #}
}

}
sub mean {
my $count=0;
my(@numbers) =@_;
foreach (@_) { $count+=$_;}
return $count/(scalar @_);
}

sub total {
my $count=0;
my(@numbers) =@_;
foreach (@_) { $count+=$_;}
return $count;
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

