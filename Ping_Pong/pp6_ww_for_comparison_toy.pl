#!/usr/bin/perl
BEGIN { unshift @INC,"/home/wangw1/bin/";}
require "restrict_digts.pm";
use File::Basename;
#print "pp";
#print map {"\t$_"} @ARGV;
#print "\n";

for ($i=0; $i<$ARGV[2]; $i++) {
   #print $ARGV[$i];

   for ($j=0; $j<$ARGV[2]; $j++) {
   $file1=fileparse($ARGV[$i]);  #$file1 is the target strand
   @namefield=split(/\./,$file1);
   #$name1=$namefield[2]."_".$namefield[1];
   $name1=$namefield[2];
   $file1=$name1;
    
   $file2=fileparse($ARGV[$j]); #file2 is the guide strand
   @namefield=split(/\./,$file2);
   $name2=$namefield[2]."_".$namefield[1];
   $file2=$name2;
    #print "$file1-$file2";
    #print "Z-score\tPercentage\tMean\tStd\tP(10)\tZ-score\tPercentage\tMean\tStd\tP(10)\n";
    open PP, ">$file1_$file2.pp";
    open OUT, ">$file1_$file2.zscore.out";
    open PPSEQ, ">$file1_$file2.ppseq";
     $X=0; $Z=0; %score=();

     foreach ($n=1;$n<=32;$n++) {
     %pp=(); %pos=(); %pp_seq=(); %pos_seq=();
     open IN, $ARGV[$i];
     while(<IN>) { chomp; s/\s+/\t/g;split(/\t/);
     next if (length($_[0])>29 || length($_[0])<23);
     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
     if ($4 eq '+') { $start=$2+$n-1; $pp{"$1:$start-"}+=$_[1]/$_[6]; $pp_seq{"$1:$start-"}=$_[0] if($n==10);}
     else { $start=$3-$n+1;$pp{"$1:$start+"}+=$_[1]/$_[6];$pp_seq{"$1:$start+"}=$_[0] if($n==10);}
     }
   
     open IN, $ARGV[$j];
     while(<IN>) { chomp; s/\s+/\t/g;split(/\t/);
     next if (length($_[0])>29 || length($_[0])<23);
     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
     if ($4 eq '+') { $pos{"$1:$2+"}+=$_[1]/$_[6]; $pos_seq{"$1:$2+"}= $_[0]; }
     else { $pos{"$1:$3-"}+=$_[1]/$_[6]; $pos_seq{"$1:$3-"}=$_[0]; }
     } 
     foreach (keys %pos)
     {
         if ($_ && exists $pp{$_})
         {
         $score{$n}+=$pos{$_}*$pp{$_} ;
         print PPSEQ "$pp_seq{$_}\t$pp{$_}\t$pos_seq{$_}\t$pos{$_}\n" if ($n==10);
         }
      }
     
     $score{$n}=0 if (!exists $score{$n});
     print PP "$n\t$score{$n}\n";
     #if ($n==10) { $X=$score{$n}; delete $score{$n};}
     }
   $X=$score{10};
   $p9=$score{9};
   $p11=$score{11};
   
   #%score_backup=%score;
   #print "peak_space\tp10\tp9\tp11";
   #for($w=20;$w<=32;$w++)
   #{
   # 
   # print "\tz_score_win$w\tpercentage_win$w\tmean_win$w\tstd_win$w";
   #}
   #print "\n";
   #print "$file1\t$X\t$p9\t$p11";
   #print OUT "$X\t$p9\t$p11";
   
   for($w=20;$w<=32;$w++)
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
   print "$file1\t$X\t$Z\t$percentage\t$m\t$std\twinsize_$w\n";
   print OUT "$file1\t$X\t$Z\t$percentage\t$m\t$std\twinsize_$w\n";
   }
   }
   #print "\n";
   #print OUT "\n";
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

