#!/usr/bin/perl
use File::Basename;
use Compress::Zlib;
#print "pp";
#print map {"\t$_"} @ARGV;
#print "\n";

for ($i=0; $i<$ARGV[2]; $i++) {
   #print $ARGV[$i];

   for ($j=0; $j<$ARGV[2]; $j++) {
   $file1=fileparse($ARGV[$i]);  #$file1 is the target strand
   $file1 =~ /(.*)\.xkxh/;
   $filename1=$1;
    
   $file2=fileparse($ARGV[$j]); #file2 is the guide strand
   $file2 =~ /(.*)\.xkxh/;
   $filename2=$1;
   if($ARGV[4])
   {
   $Transposon=$ARGV[4];
   print "$Transposon\t$filename1-$filename2";
   open PPSEQ, ">$ARGV[3]/$filename1_$filename2.$Transposon.ppseq";
   }
   else
   {print "$filename1-$filename2";
   	open PPSEQ, ">$ARGV[3]/$filename1_$filename2.ppseq";
   }
    
    
     $X=0; $Z=0; %score=();

     foreach ($n=1;$n<=20;$n++) {
     %pp=(); %pos=(); %pp_seq=(); %pos_seq=();
     if($filename1=~/gz/)
     {
	     $gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
	     while($gz->gzreadline($_) > 0)
	          { chomp; split(/\t/);
		     next if (length($_[0])>29 || length($_[0])<23);
		     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
		     if ($4 eq '+') { $start=$2+$n-1; $pp{"$1:$start-"}+=$_[1]/$_[6]; $pp_seq{"$1:$start-"}=$_[0];}
		     else { $start=$3-$n+1;$pp{"$1:$start+"}+=$_[1]/$_[6];$pp_seq{"$1:$start+"}=$_[0];}
		     }
     }
     else
     {
	     open IN, $ARGV[$i];
	     while(<IN>) 
	          { chomp; split(/\t/);
			     next if (length($_[0])>29 || length($_[0])<23);
			     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
			     if ($4 eq '+') { $start=$2+$n-1; $pp{"$1:$start-"}+=$_[1]/$_[6]; $pp_seq{"$1:$start-"}=$_[0];}
			     else { $start=$3-$n+1;$pp{"$1:$start+"}+=$_[1]/$_[6];$pp_seq{"$1:$start+"}=$_[0];}
			  }
     }

   
      if($filename1=~/gz/)
     {
	     $gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
	     while($gz->gzreadline($_) > 0)
	          { chomp; split(/\t/);
			     next if (length($_[0])>29 || length($_[0])<23);
			     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
			     if ($4 eq '+') { $pos{"$1:$2+"}+=$_[1]/$_[6]; $pos_seq{"$1:$2+"}= $_[0]; }
			     else { $pos{"$1:$3-"}+=$_[1]/$_[6]; $pos_seq{"$1:$3-"}=$_[0]; }
			     } 
     }
     else
     {
	     open IN, $ARGV[$i];
	     while(<IN>) 
	          { chomp; split(/\t/);
		     next if (length($_[0])>29 || length($_[0])<23);
		     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
		     if ($4 eq '+') { $pos{"$1:$2+"}+=$_[1]/$_[6]; $pos_seq{"$1:$2+"}= $_[0]; }
		     else { $pos{"$1:$3-"}+=$_[1]/$_[6]; $pos_seq{"$1:$3-"}=$_[0]; }
		     } 
     }

     foreach (keys %pos)
     {
         if ($_ && exists $pp{$_})
         {
         $score{$n}+=$pos{$_}*$pp{$_} ;
         print PPSEQ "$pp_seq{$_}\t$pp{$_}\t$pos_seq{$_}\t$pos{$_}\n" if ($n==10);
         }
      }
     #print "$n\t$score{$n}\n";
     $score{$n}=0 if (!exists $score{$n});
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