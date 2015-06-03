#!perl -slw
use strict;
use threads;
use threads::shared;
use Thread::Queue;
use Data::Dump qw[ pp ];

sub helper {
    my $Q = shift;
    while( my $ref = $Q->dequeue ) {;
        lock $ref;
        $ref->{NEW_KEY} = 1;
    }
}

sub my_sub {
    my( $ref, $n ) = @_;
    my $Q = new Thread::Queue;
    my @threads = map async( \&helper, $Q ), 1 .. $n;
    $Q->enqueue( values %{ $ref } );
    $Q->enqueue( (undef) x $n );
    $_->join for @threads;
}

my $hoh = {
    A => shared_clone( { NAME => 'aa' } ),
    B => shared_clone( { NAME => 'bb' } ),
};

pp $hoh;
my_sub( $hoh, 2 );
pp $hoh;