#!/usr/bin/perl

while(<>) { chomp; split(/\t/);
$_[2]=~/(\w.+):(\d+)-(\d+)\((.+)\)/;
$key="$1\t$2\t$3\t$4\t$_[0]\t$_[1]\t$_[6]\t$_[3]\n";
$hash{$key}=0;
}

foreach (keys %hash) { 
print ; 
}
