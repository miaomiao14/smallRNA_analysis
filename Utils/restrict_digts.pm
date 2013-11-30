#$a=&restrict_num_decimal_digits(3.1245,3);
#print "$a\n";

sub restrict_num_decimal_digits
{
  my $num=shift;#the number to work on
  my $digs_to_cut=shift;# the number of digits after
		  	    # the decimal point to cut
		#(eg: $digs_to_cut=3 will leave
		# two digits after the decimal point)

  if ($num=~/\d+\.(\d){$digs_to_cut,}/)
  {
    # there are $digs_to_cut or
    # more digits after the decimal point
    $num=sprintf("%.".($digs_to_cut-1)."f", $num);
  }
  return $num;
}
1;