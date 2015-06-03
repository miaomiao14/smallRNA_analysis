#!/usr/bin/perl
use strict;
use warnings;
# array reclaims mem 

show_size();

{
    my @var=(0..1000000);
    show_size();
    undef @var;
}

show_size();
exit(0);

sub show_size {
    local $/;
    open(my $pfh, '<', "/proc/$$/status") || die $!;
    my $size = <$pfh> =~ /VmSize:\s+(\d+)/
               ? $1
               : 'unknown';
    close($pfh);
    print "Process size: $size\n";
}

