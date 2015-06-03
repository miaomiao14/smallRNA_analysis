#!perl -slw
use strict;
use threads;
use threads::shared;
use Data::Rmap;
use Data::Dumper;
$Data::Dumper::Terse=1; 
$Data::Dumper::Indent=0;

sub thread {
    my( $hashref ) = @_;
    rmap{
        print "$_";
		$_ = uc $_ ; 
    } $hashref;
}

my %hash : shared;

$hash{ leaf } = 'fred';
$hash{ L1 } = &share( {} );
$hash{ L1 }{ leaf } = 'wilma';
$hash{ L1 }{ L2 } = &share( {} );
$hash{ L1 }{ L2 }{ leaf } = 'bam bam';
$hash{ L1 }{ L2 }{ L3 } = &share( {} );
$hash{ L1 }{ L2 }{ L3 }{ leaf } = 'flintstone';

threads->create( \&thread, \%hash )->join;
print Dumper(%hash);

__END__


# You can only assign references to shared arrays and hashes into shared hashes. You can either declare them shared and assign a reference:
#
# my %hash : shared;
# my @array : shared;
# my %hash2 : shared;
#
# $hash{ array } = \@array;
# $hash{ hash  } = \%hash;
# [download]
# Which is simple but often inconvenient. Or you can use threads::shared::share to them before assignment.
#
# There are two caveats with the second approach.
#
# If you share a pre-existing hash or array that already contains data, that data is discarded.
# The share() sub uses a prototype which rejects attempts to share anonymous structures.
# You can make life a little easier by bypassing the prototypes using $sharedHash{ key } = &share( {} ) and $sharedHash{ key } = &share( [] ), but don't be tempted to do $sharedHash{ key } = &share( [ 'some', 'stuff', 'here' ] );because the contents will be discarded.