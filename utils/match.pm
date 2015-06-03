#!/usr/bin/perl
sub match {
	
	my(@files) =@_;
	open IN, $files[0];
	while(my $line=<IN>) {
	chomp $line;
	my @l=split(/\t/,$line);
	$hash{$l[0]}=0;
	}

	open IN, $files[1];
	while(my $line=<IN>) {
	chomp $line;
	@l=split(/\t/,$line);
	print "$line" if (exists $hash{$l[0]});
	}
}
