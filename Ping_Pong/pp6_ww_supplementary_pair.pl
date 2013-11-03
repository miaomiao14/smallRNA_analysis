#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
#print "pp";
#print map {"\t$_"} @ARGV;
#print "\n";

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
    open PPSEQ, ">$file1_$file2.ppseq";
     $X=0; $Z=0; %score=();%s=();

     foreach ($n=-50;$n<=50;$n++) {
     %pp=(); %pos=(); %pp_seq=(); %pos_seq=();
     open IN, $ARGV[$i];
     while(<IN>) { chomp; split(/\t/);
     next if (length($_[0])>29 || length($_[0])<23);
     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
     if ($4 eq '+') { $start=$2+$n-1; $pp{"$1:$start-"}+=$_[1]/$_[6]; $pp_seq{"$1:$start-"}=$_[0];}
     else { $start=$3-$n+1;$pp{"$1:$start+"}+=$_[1]/$_[6];$pp_seq{"$1:$start+"}=$_[0];}
     }
   
     open IN, $ARGV[$j];
     while(<IN>) { chomp; split(/\t/);
     next if (length($_[0])>29 || length($_[0])<23);
     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
     if ($4 eq '+') { $pos{"$1:$3+"}+=$_[1]/$_[6]; $pos_seq{"$1:$3+"}= $_[0]; }  ##use end instead of start
     else { $pos{"$1:$2-"}+=$_[1]/$_[6]; $pos_seq{"$1:$2-"}=$_[0]; }  ## use end instead of start
     } 
     foreach (keys %pos)
     {
         if ($_ && exists $pp{$_})
         {
         $score{$n}+=$pos{$_}*$pp{$_} ;
         $s{$n}+=&min($pos{$_},$pp{$_});
         $spaces=length($pos_seq{$_})-10;
         #printf PPSEQ "%$spaces";
         print PPSEQ "$pp_seq{$_}\t$pp{$_}\n$pos_seq{$_}\t$pos{$_}\n" if ($n==10);
         }
      }
     #print "$n\t$score{$n}\n";
     $score{$n}=0 if (!exists $score{$n});
     print PP1 "$n\t$score{$n}\n";
     print PP2 "$n\t$s{$n}\n";
     if ($n==10) { $X=$score{$n}; delete $score{$n};}
     }
   $std=&standard_deviation(values %score);
   if ($std>0) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
   print "\t$Z\t$X\n";
   }
   #print "\n";
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

