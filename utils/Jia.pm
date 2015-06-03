
my %RevCompBasePairs=qw/A T T A G C C G a t t a g c c g U A u a R Y r y Y R y r M K m k K M k m S S s s W W w w H D D H h d d h B V V B b v v b N N n n/;
sub revfa {
  my $s=shift @_;
  my $new_s='';
  for my $b (split//,$s) {
   if(!exists $RevCompBasePairs{$b}) {$new_s="?".$new_s;}
   else {$new_s=$RevCompBasePairs{$b}.$new_s;}
  }
  return $new_s;
}

sub permute {
my $seq1=$_[0];
my $seq2=$_[1];
my $i;
if (length($seq2)<=1) { print $seq1.$seq2,"\n";}
else {
    for ($i=0; $i<length($seq2);$i++) {
    @seq2=split(//,$seq2);
    $new_seq2=join('',@seq2[0..$i-1,$i+1..length($seq2)-1]);
    $new_seq1=$seq1.$seq2[$i];
    &permute($new_seq1,$new_seq2);
    }
}
}


1;
