#$transPairSuppSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}{$diff->list}+=1/$NTM{$l[2]};
#$transPairSuppReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}{$diff->list}+=$guideQueryReads*$targetpiGuideReads/$NTM{$l[2]};
use PDL;
use PDL::Char;
#$PDL::SHARE = $PDL::SHARE; # keep stray warning quiet
use Data::Dumper;
use strict;
use warnings;

#my %transPairSuppSpeciesTotal=(PDL::Char->new("011110") => 1, PDL::Char->new("011111") => 2);
#my %transPairSuppReadsTotal=(PDL::Char->new("011110") => 2.5, PDL::Char->new("011111") => 10);
my $ref=("CAAAACT");
my $source = PDL::Char->new($ref);
my @query=("AAAAACT","AAAAGCT");
my %transPairSuppSpeciesTotal;
my %transPairSuppReadsTotal;
my @diffs;
my $i=0;
for my $str (@query) {
	my $match = PDL::Char->new($str);
	my $diff = $match == $source;
	#push @diffs, $diff->list;
	my $difftemp=join('',$diff->list);
	#print $difftemp,"\n";
	$transPairSuppSpeciesTotal{ "$difftemp" } =1.5;
	$transPairSuppReadsTotal{ "$difftemp" } =2;	
	$i++;
	}


#my $ref1=PDL::Char->new("011110");
#my $ref2=PDL::Char->new("011110");
#my %transPairSuppSpeciesTotal=($diffs[0] => 1, $diffs[1] => 2);
#my %transPairSuppReadsTotal=($diffs[0]  => 2.5, $diffs[1] => 10);

my ($scaledSpeciesRef,$scaledReadsRef)=&SuppComBitSum(\%transPairSuppSpeciesTotal,\%transPairSuppReadsTotal);

for(my $position=0; $position< @{$scaledSpeciesRef};$position++)
{
	
	print "Speciesindex:$position\t$scaledSpeciesRef->[$position]\n";
}
for(my $position=0; $position< @{$scaledReadsRef};$position++)
{
	
	print "Readsindex:$position\t$scaledReadsRef->[$position]\n";
}


sub SuppComBitSum
{
	my ($transPairSuppSpeciesTotalRef,$transPairSuppReadsTotalRef)=@_;
	my @scaledSuppSpeciesCom=();
	my @scaledSuppReadsCom=();
	my %transPairSuppSpeciesTotal=%{$transPairSuppSpeciesTotalRef};
	my %transPairSuppReadsTotal=%{$transPairSuppReadsTotalRef};
	foreach my $bitValue (keys %transPairSuppSpeciesTotal)
	{
		print "bits:",$bitValue,"\n";
		#print Dumper $bitValue,"\n";
		my $scaleFactorSpecies=$transPairSuppSpeciesTotal{$bitValue};
		my $scaleFactorReads=$transPairSuppReadsTotal{$bitValue};
		my $i=0;
		#my @bits=split(//,unpack("b*",$bitValue));
		my @bits=split(//,$bitValue);
		print "scaleFactor:\t",$scaleFactorSpecies,"\t",$scaleFactorReads,"\n";
		foreach my $b ( @bits ) # $bitValue is a list
		{
		
			print $i,"\t",$b,"\n";
			$scaledSuppSpeciesCom[$i]+=$b*$scaleFactorSpecies;
			$scaledSuppReadsCom[$i]+=$b*$scaleFactorReads;
			
			print $scaledSuppSpeciesCom[$i],"\t",$scaledSuppReadsCom[$i],"\n";
			$i++;
		}

	}

	return \@scaledSuppSpeciesCom,\@scaledSuppReadsCom;
}




