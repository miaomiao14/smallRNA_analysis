#!/usr/bin/perl

BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
BEGIN { unshift @INC,"/home/wangw1/bin/";}
require "sort_hash_key.pm";
require "restrict_digts.pm";
use File::Basename;
use Compress::Zlib;

for ($i=0; $i<$ARGV[2]; $i++) {
   #print $ARGV[$i];

   for ($j=0; $j<$ARGV[2]; $j++) {
   $file1=fileparse($ARGV[$i]);  #$file1 is the target strand
   @namefield=split(/\./,$file1);
   $name1=$namefield[0]."_".$namefield[2]."_".$namefield[1];
   $file1=$name1;
    
   $file2=fileparse($ARGV[$j]); #file2 is the guide strand
   @namefield=split(/\./,$file2);
   $name2=$namefield[0]."_".$namefield[2]."_".$namefield[1];
   $file2=$name2;
    #print "$file1-$file2";
    #print "Z-score\tPercentage\tMean\tStd\tP(10)\tZ-score\tPercentage\tMean\tStd\tP(10)\n";
    open PP, ">$file1_$file2.pp";
    open OUT, ">$file1_$file2.zscore.out";
    #open PPSEQ, ">$file1_$file2.ppseq";
     $X=0; $Z=0; %score=();
      %pp_pos=();
     foreach ($n=1;$n<=50;$n++) {
     %pp=(); %pos=(); %pp_seq=(); %pos_seq=();
     my $gz = gzopen($ARGV[$i], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
     #open IN, $ARGV[$i];
     while($gz->gzreadline($_) > 0) { chomp; s/\s+/\t/g;split(/\t/);
     next if (length($_[0])>29 || length($_[0])<23);
     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
     if ($4 eq '+') { $start=$2+$n-1; $pp{"$1:$start-"}+=$_[1]/$_[6]; }
     else { $start=$3-$n+1;$pp{"$1:$start+"}+=$_[1]/$_[6];}
     }
      #die "Error reading from $ARGV[$i]: $gzerrno\n" if $gzerrno != Z_STREAM_END ;
      $gz->gzclose() ;
      
     my $gz2 = gzopen($ARGV[$j], "rb") or die "Cannot open $ARGV[$j]: $gzerrno\n" ;
     while($gz2->gzreadline($_) > 0) { chomp; s/\s+/\t/g;split(/\t/);
     next if (length($_[0])>29 || length($_[0])<23);
     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
     if ($4 eq '+') { $pos{"$1:$2+"}+=$_[1]/$_[6];  }
     else { $pos{"$1:$3-"}+=$_[1]/$_[6];  }
     }
     #die "Error reading from $ARGV[$j]: $gzerrno\n" if $gzerrno != Z_STREAM_END ;
     $gz2->gzclose() ;
     foreach (keys %pos)
     {
         if ($_ && exists $pp{$_})
         {
         $score{$n}+=$pos{$_}*$pp{$_} ;
         $_=~/(chr.+):(\d+)(\+)/;
         $pp_pos{$1}{$3}{$2}=1 if($n==10 && $1);
         }
      }
     
     $score{$n}=0 if (!exists $score{$n});
     print PP "$n\t$score{$n}\n";
     }
   %crowd=();
   $total_pos=0;
   foreach $chr (keys %pp_pos)
   {
      foreach $strand (keys %{$pp_pos{$chr}})
      {
         %ppchrpos=%{$pp_pos{$chr}{$strand}};
         #@ppchrpos_sort=&sort_hash_key(%ppchrpos);
         $total_pos+=scalar( keys %ppchrpos);
         $last_pos=0;
         foreach $p (sort { $a <=> $b } keys %ppchrpos)
         {
            if($p>$last_pos)
            {
               for($i=1;$i<10000000;$i++)
               {
                  $current_pos=$p+$i;
                  if($ppchrpos{$current_pos})
                  {
                     $k=0;
                  }
                  else
                  {
                     #if($i>1)
                     #{
                        $crowdedness=$i;
                        $crowd{$crowdedness}++;
                     #}
		     $last_pos=$p+$i-1;
                     last;
                  }
               }
            }
         }
         #last; ##no need to calculate twice for + and -
         
      }
   }
   open CROWD, ">$file1_$file2.crowdedness.out";
   foreach $crowd (sort { $a <=> $b } keys %crowd)
   {
      print CROWD "$crowd\t$crowd{$crowd}\n";
   }
   print CROWD "0\t$total_pos\n";
   close(CROWD);
   $X=$score{10};
   $p9=$score{9};
   $p11=$score{11};
   for($w=11;$w<=50;$w++)
   {
   #$hash_name="s_".$w;
   %s=();
   for($j=1;$j<=$w;$j++)
   {
      $s{$j}=$score{$j};
   }
  
   $total_score=0;
   $total_score=&total(values %s);
   if ($total_score!=0)
   {
   $percentage=$X/$total_score;
   $percentage=&restrict_num_decimal_digits($percentage,5);
   }
   else
   {
      $percentage=0;
   }
   delete $s{10};
   $std=&standard_deviation(values %s);
   $std=&restrict_num_decimal_digits($std,3);
   $m=&mean(values %s);
   $m=&restrict_num_decimal_digits($m,3);
   
   if ($X!=0 && $std>0) { $Z=($X-&mean(values %s))/$std; $Z=&restrict_num_decimal_digits($Z,4);} else {$Z=-10;}
   #print "$file1\t$X\t$Z\t$percentage\t$m\t$std\twinsize_$w\n";
   print OUT "$file1\t$X\t$Z\t$percentage\t$m\t$std\twinsize_$w\n";
   }
   
   }
   #print "\n";
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

