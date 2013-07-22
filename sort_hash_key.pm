#%test=("10"=>1,"11"=>2,"9"=>3);

#&sort_hash_key(%test);
sub sort_hash_key
{
my %hash=@_;
my @sortarray=();
foreach $k (sort {$a <=> $b } keys %hash)
{
  push @sortarray, $k;
}
#print "@sortarray";
return @sortarray;	
}

#%t1=("110"=> 10, "102"=> 11);
#$t1ref=\%t1;
#%t2=("210"=> 10, "202"=> 11);
#$t2ref=\%t2;
#%test=("280" => $t1ref,"80" => $t2ref);
#foreach $start ( keys %test)
#{
#	foreach $end (keys %{$test{$start}})
#	{
#		print "$start\t$end\t$test{$start}{$end}\n";
#	}
#}
#&sort_double_key_hash(%test);

sub sort_double_key_hash
{
my %hash= @_;
my @sortarray=();
foreach $start ( sort {$a <=> $b } keys %hash)
{
    foreach $end ( sort {$a <=> $b} keys %{$hash{$start}})
    {
        #print "$start\t$end\t$hash{$start}{$end}\n";
        push @sortarray, [$start,$end,$hash{$start}{$end}];
    }
}
return @sortarray;
}

#@test=([1,2,3],[241,2,3],[31,2,3]);
#&sort_arrayref_value(@test);
#@test=([1,2,3]);
#&sort_arrayref_value(@test);
sub sort_arrayref_value
{
my @array= @_;
my @sortarray=();
foreach $k ( sort {$a->[0] <=> $b->[0] } @array)
{
#	print "$k->[0],$k->[1],$k->[2]\n";
	push @sortarray, [$k->[0],$k->[1],$k->[2]];
}
return @sortarray;
}

sub rsort_arrayref_value
{
my @array= @_;
my @sortarray=();
foreach $k ( sort {$b->[0] <=> $a->[0] } @array)
{
        #print "$k->[0],$k->[1],$k->[2]\n";
        push @sortarray, [$k->[0],$k->[1],$k->[2]];
}
return @sortarray;
}

sub sort_arrayref2_value
{
my @array= @_;
my @sortarray=();
foreach $k ( sort {$a->[0] <=> $b->[0] } @array)
{
        #print "$k->[0],$k->[1],$k->[2]\n";
        push @sortarray, [$k->[0],$k->[1]];
}
return @sortarray;
}

1;
