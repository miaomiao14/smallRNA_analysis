#!/usr/bin/perl
use strict;
use warnings;
use MLDBM qw(DB_File);
use Fcntl;

my %hash;
my %inner;
my $id=$ARGV[0];
my $file=$ARGV[1];
my $counter=1;
tie (%hash,"MLDBM","$id.db",O_CREAT|O_RDWR, 0640) or die "$!\n";

open FILE, "<$file" or die "$!\n";
while(<FILE>){
	if($. == 1){
		next;
	}
	else{
		chomp;
		my @line=split(/,/,$_);
		for(my $i=0;$i<scalar(@line);$i++){ 
			$hash{$counter}=$line[$i];
			$counter++;
		}
	}
}
close FILE;

