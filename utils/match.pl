#!/usr/bin/perl

open IN, $ARGV[0];
while(my $line=<IN>) {
chomp $line;
my @l=split(/\t/,$line);
$hash{$l[0]}=0;
}

open IN, $ARGV[1];
while(my $line=<IN>) {
chomp $line;
@l=split(/\t/,$line);
print "$line\n" if (exists $hash{$l[0]});
}
