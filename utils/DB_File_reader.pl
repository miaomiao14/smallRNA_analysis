#!/usr/bin/perl
use strict;
use warnings;
use MLDBM qw(DB_File);
use Fcntl;

my %hash;
my %inner;
my $counter=1;
tie (%hash,"MLDBM","DNA_shape_db/rv_helt.db") or die "$!\n";



for(my $i=1;$i<1001;$i++){
	print "$hash{$i}\n";
}


