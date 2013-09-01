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
   open PPSEQ, ">$ARGV[3]/$filename1.$filename2.$Transposon.ppseq";
   open PPSCORE, ">$ARGV[3]/$filename1.$filename2.$Transposon.ppscore";
   #print PPSCORE, "$Transposon\t$filename1-$filename2";
   }
   else
   {print "$filename1-$filename2";
   	open PPSEQ, ">$ARGV[3]/$filename1.$filename2.ppseq";
   	open PPSCORE, ">$ARGV[3]/$filename1.$filename2.ppscore";
   	#print PPSCORE "$filename1-$filename2";
   }
   #total transposon piRNAs
   $totalReads1=0;
   $totalReads2=0;
   $ppRead1=0;
   $ppRead2=0;
   if($file1=~/gz/)
   {
    $gz = gzopen($ARGV[$i], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
	     while($gz->gzreadline($_) > 0)
	     {chomp; split(/\t/); $totalReads1+=$_[1]/$_[6]; }
	 $gz->gzclose();  
   }
   else
   {
   	 open IN, $ARGV[$i] or die "Cannot open $ARGV[$i]: $!\n";
	     while(<IN>) 
	     { chomp; split(/\t/);$totalReads1+=$_[1]/$_[6]; }
	     close(IN);
   	
   }
   
   if($file2=~/gz/)
   {
    $gz = gzopen($ARGV[$j], "rb") or die "Cannot open $ARGV[$j]: $gzerrno\n" ;
	     while($gz->gzreadline($_) > 0)
	     {chomp; split(/\t/); $totalReads2+=$_[1]/$_[6]; }
	 $gz->gzclose();  
   }
   else
   {
   	 open IN, $ARGV[$j] or die "Cannot open $ARGV[$j]: $!\n";
	     while(<IN>) 
	     { chomp; split(/\t/);$totalReads2 +=$_[1]/$_[6]; }
	     close(IN);
   	
   }
   
   
   	
    
     $X=0; $Z=0; %score=();

     foreach ($n=1;$n<=20;$n++) {
     %pp=(); %pos=(); %pp_seq=(); %pos_seq=();
     if($file1=~/gz/)
     {
	     $gz = gzopen($ARGV[$i], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
	     while($gz->gzreadline($_) > 0)
	          { chomp; split(/\t/);
	         	   
		     #next if (length($_[0])>29 || length($_[0])<23);
		     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
		     if ($4 eq '+') { $start=$2+$n-1; $pp{"$1:$start-"}+=$_[1]/$_[6]; $pp_seq{"$1:$start-"}=$_[0];}
		     else { $start=$3-$n+1;$pp{"$1:$start+"}+=$_[1]/$_[6];$pp_seq{"$1:$start+"}=$_[0];}
		     }
		  $gz->gzclose();
     }
     else
     {
	     open IN, $ARGV[$i] or die "Cannot open $ARGV[$i]: $!\n";
	     while(<IN>) 
	          { chomp; split(/\t/);
			     #next if (length($_[0])>29 || length($_[0])<23);
			     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/;
			     if ($4 eq '+') { $start=$2+$n-1; $pp{"$1:$start-"}+=$_[1]/$_[6]; $pp_seq{"$1:$start-"}=$_[0];}
			     else { $start=$3-$n+1;$pp{"$1:$start+"}+=$_[1]/$_[6];$pp_seq{"$1:$start+"}=$_[0];}
			  }
		 close(IN);
     }

   
      if($file2=~/gz/)
     {
	     $gz = gzopen($ARGV[$j], "rb") or die "Cannot open $ARGV[$j]: $gzerrno\n" ;
	     while($gz->gzreadline($_) > 0)
	          { chomp; split(/\t/);
			     #next if (length($_[0])>29 || length($_[0])<23);
			     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
			     if ($4 eq '+') { $pos{"$1:$2+"}+=$_[1]/$_[6]; $pos_seq{"$1:$2+"}= $_[0]; }
			     else { $pos{"$1:$3-"}+=$_[1]/$_[6]; $pos_seq{"$1:$3-"}=$_[0]; }
			     } 
		  $gz->gzclose();
     }
     else
     {
	     open IN, $ARGV[$j]  or die "Cannot open $ARGV[$j]: $!\n";
	     while(<IN>) 
	          { chomp; split(/\t/);
		     #next if (length($_[0])>29 || length($_[0])<23);
		     $_[2]=~/(chr.+):(\d+)-(\d+)\((.+)\)/; 
		     if ($4 eq '+') { $pos{"$1:$2+"}+=$_[1]/$_[6]; $pos_seq{"$1:$2+"}= $_[0]; }
		     else { $pos{"$1:$3-"}+=$_[1]/$_[6]; $pos_seq{"$1:$3-"}=$_[0]; }
		     } 
		  close(IN);
     }

     foreach (keys %pos)
     {
         if ($_ && exists $pp{$_})
         {
         $score{$n}+=$pos{$_}*$pp{$_} ;
         print PPSEQ "$pp_seq{$_}\t$pp{$_}\t$pos_seq{$_}\t$pos{$_}\n" if ($n==10);
         $ppRead1+=$pp{$_} if ($n==10);
         $ppRead2+=$pos{$_} if ($n==10);
         }
      }
     print PPSCORE "$n\t$score{$n}\n";
     $score{$n}=0 if (!exists $score{$n});
     if ($n==10) { $X=$score{$n}; delete $score{$n};}
     }
   $std=&standard_deviation(values %score);
   if ((scalar %score)>9 && $std>0) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
   $ppRatio1=$ppRead1/$totalReads1*100;
   $ppRatio2=$ppRead2/$totalReads2*100;
   $roundedppRatio1 = sprintf("%.2f", $ppRatio1);
   $roundedppRatio2 = sprintf("%.2f", $ppRatio2);
   $roundedZ = sprintf("%.2f", $Z);
   $roundedX = sprintf("%.2f", $X);
   $roundedppReads1 = sprintf("%.2f", $ppRead1);
   $roundedppReads2 = sprintf("%.2f", $ppRead2);
   $roundedtotalReads1 = sprintf("%.2f", $totalReads1);
   $roundedtotalReads2 = sprintf("%.2f", $totalReads2);
   print "\t$roundedZ\t$roundedX\t$roundedppReads1\t$roundedtotalReads1\t$roundedppRatio1%\t$roundedppReads2\t$roundedtotalReads2\t$roundedppRatio2%\n";

   }
   #print "\n";
}
close(PPSCORE);
close(PPSEQ);

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