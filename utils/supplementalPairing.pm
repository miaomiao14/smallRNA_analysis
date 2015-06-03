sub SuppComBitSum
{
	my ($transPairSuppSpeciesTotalRef,$transPairSuppReadsTotalRef)=@_;
	my @scaledSuppSpeciesCom=();
	my @scaledSuppReadsCom=();
	my %transPairSuppSpeciesTotal=%{$transPairSuppSpeciesTotalRef};
	my %transPairSuppReadsTotal=%{$transPairSuppReadsTotalRef};
	foreach my $bitValue (keys %transPairSuppSpeciesTotal)
	{
		#print "bits:",$bitValue,"\n";
		#print Dumper $bitValue,"\n";
		my $scaleFactorSpecies=$transPairSuppSpeciesTotal{$bitValue};
		my $scaleFactorReads=$transPairSuppReadsTotal{$bitValue};
		my $i=1; #starting postion should be seed +1;
		#my @bits=split(//,unpack("b*",$bitValue));
		my @bits=split(//,$bitValue);
		#print "scaleFactor:\t",$scaleFactorSpecies,"\t",$scaleFactorReads,"\n";
		foreach my $b ( @bits ) # $bitValue is a list
		{

			#print $i,"\t",$b,"\n";
			if($b==0)
			{
				last; #if meet any 0, then it's not continuous; by doing this, will miss a lot of binary string
			}
			else
			{
				$scaledSuppSpeciesCom[$i]+=$scaleFactorSpecies; #position's index increase and add 1 
				$scaledSuppReadsCom[$i]+=$scaleFactorReads;
				$i++;
				#print $scaledSuppSpeciesCom[$i],"\t",$scaledSuppReadsCom[$i],"\n";
			}
		}

	}

	return \@scaledSuppSpeciesCom,\@scaledSuppReadsCom;

}


sub SuppComContinuousBitSum
{
	my ($transPairSuppSpeciesTotalRef,$transPairSuppReadsTotalRef)=@_;
	my @scaledSuppSpeciesCom=();
	my @scaledSuppReadsCom=();
	my %transPairSuppSpeciesTotal=%{$transPairSuppSpeciesTotalRef};
	my %transPairSuppReadsTotal=%{$transPairSuppReadsTotalRef};
	foreach my $bitValue (keys %transPairSuppSpeciesTotal)
	{
		#print "bits:",$bitValue,"\n";
		#print Dumper $bitValue,"\n";
		my $scaleFactorSpecies=$transPairSuppSpeciesTotal{$bitValue};
		my $scaleFactorReads=$transPairSuppReadsTotal{$bitValue};
		my $i=0;
		#my @bits=split(//,unpack("b*",$bitValue));
		my @bits=split(//,$bitValue);
		#print "scaleFactor:\t",$scaleFactorSpecies,"\t",$scaleFactorReads,"\n";
		foreach my $b ( @bits ) # $bitValue is a list
		{

			#print $i,"\t",$b,"\n";
			$scaledSuppSpeciesCom[$i]+=$b*$scaleFactorSpecies;
			$scaledSuppReadsCom[$i]+=$b*$scaleFactorReads;
			$i++;
			#print $scaledSuppSpeciesCom[$i],"\t",$scaledSuppReadsCom[$i],"\n";
		}

	}

	return \@scaledSuppSpeciesCom,\@scaledSuppReadsCom;

}




1;