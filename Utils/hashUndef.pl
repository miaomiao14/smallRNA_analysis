show_size();

{
    my %hash;

    for(0..1000000){
      $hash{$_}= 'perlperlperl';
    
    }
    
    show_size();
    
    foreach(0..1000000){
      $hash{$_}= '';
      undef $hash{$_};
    }
    
    undef %hash;
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


