sub split_chr {
	my ($file,$OUTDIR)=@_;
	my %filehash=();
	if($file=~/gz/)
	{
		my $gz="";
		$gz = gzopen($file, "rb") or die "Cannot open $file: $gzerrno\n" ;
		while($gz->gzreadline(my $record) > 0)
		{ 
			chomp $record;
			my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
			$filehash{$chr}{$record}=1;
		}
		$gz->gzclose();
	}
	else
	{
		open IN, $file or die "Cannot open $file: $!\n";
		while(my $record=<IN>) 
		{
			chomp $record;
			my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
			$filehash{$chr}{$record}=1;
		}
		close(IN);
	}
	foreach my $chr (keys %filehash)
	{
		open OUT, ">$OUTDIR/$file.$chr" or die "Cannot open $OUTDIR/$file.$chr to write: $!\n";
		foreach my $record (keys %{$filehash{$chr}})
		{
			print OUT $record,"\n";
		}
		close(OUT);
	}
		
}

1;