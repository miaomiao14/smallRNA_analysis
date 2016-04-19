#!/usr/bin/perl

     $file1=$ARGV[0]; $file2=$ARGV[0];
     open OUT, "> $ARGV[0].pp";
     $X=0; $Z=0; %score=();

     foreach ($n=1;$n<=20;$n++) {
     %pp=(); %pos=();
     open IN, $file1;
     while(<IN>) { chomp; split(/\t/);
     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
     if ($4 eq '+') { $start=$2+$n-1; $pp{"$1:$start-"}+=$_[1]/$_[6]; }
     else { $start=$3-$n+1;$pp{"$1:$start+"}+=$_[1]/$_[6];}
     }
   
     open IN, $file2;
     while(<IN>) { chomp; split(/\t/);
     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
     if ($4 eq '+') { $pos{"$1:$2+"}+=$_[1]/$_[6];}
     else { $pos{"$1:$3-"}+=$_[1]/$_[6];}
     } 
     foreach (keys %pos) {  $score{$n}+=$pos{$_}*$pp{$_} if ($_ && exists $pp{$_});}
     $score{$n}=0 if (!exists $score{$n});
     print OUT "$n\t$score{$n}\n";
     if ($n==10) { $X=$score{$n}; delete $score{$n};}
     }
   $std=&standard_deviation(values %score);
   if ($std>0) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
   print $Z;


sub mean {
my $count=0;
my(@numbers) =@_;
foreach (@_) { $count+=$_;}
return $count/(scalar @_);
}

sub standard_deviation {
my(@numbers) = @_;
#Prevent division by 0 error in case you get junk data
return () unless(scalar(@numbers));

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

