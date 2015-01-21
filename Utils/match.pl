#!/usr/bin/perl

open IN, $ARGV[0];
while(<IN>) {
chomp;
split(/\t/);
$hash{$_[0]}=0;
}

open IN, $ARGV[1];
while(<IN>) {
chomp;
split(/\t/);
print "$_\n" if (exists $hash{$_[0]});
}
