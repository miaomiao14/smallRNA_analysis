#Jia Xu
#List of subroutines
# mean
# std
# min
# max

sub mean {
my $count=0;
my(@numbers) =@_;
foreach (@_) { $count+=$_;}
if (scalar @_ >0) {return $count/(scalar @_);} else { return 0;}
}

sub std {
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

# min(x,y) returns minimun of x and y.
sub min {
my $min=$_[0]>$_[1]?$_[1]:$_[0];
my $i;
 foreach ($i=1;$i<@_;$i++) {
 $min=$_[$i]>$min?$min:$_[$i];
}
return $min;
}

# max(x,y) returns maximun of x and y.
sub max {
my $max=$_[0]>$_[1]?$_[0]:$_[1];
my $i;
 foreach ($i=1;$i<@_;$i++) {
 $max=$_[$i]<$max?$max:$_[$i];
}
return $max;
}

1;
