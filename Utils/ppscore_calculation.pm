use strict;
use warnings;
use File::Basename;
use Compress::Zlib;


	my $file1=$ARGV[0];
	my $file2=$ARGV[1];
	my $filetemp1="";
	my $indir="";
	($filetemp1,$indir)=fileparse($file1);
	$filetemp1 =~ /(.*)\.bed*.*/;
	my $filename1="";
	$filename1=$1;

	my $filetemp2="";
	$filetemp2=fileparse($file2);
	$filetemp2 =~ /(.*)\.bed*.*/;
	my $filename2="";
	$filename2=$1;
	
	
	my $format=$ARGV[2];
	
	
	my $chromosome="";
	if($ARGV[3])
	{
		$chromosome=$ARGV[3];
		open PPSCORE, ">$indir/$filename1.$filename2.$chromosome.ppscore";
	}
	else
	{
		open PPSCORE, ">$indir/$filename1.$filename2.ppscore";
	}
	
	my %score=();
	
	

	
	
	foreach (my $n=1;$n<=20;$n++) #iterate with different 5'-5' distance(overlap) of piRNAs from oposite strand
	{
		my %pp=(); my %pos=(); my %pp_seq=();my %pos_seq=();
		my %pphash=(); my %poshash=();


		#the target strand	
		open IN, $file1 or die "Cannot open $file1: $!\n";
		while(my $record=<IN>) 
		{
			chomp $record;
			my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
			my $len=0;
			$len=$end-$start;
			next if ($len>29 || $len<23);
			my $ppStart=0;

			if ($strand eq '+') 
			{
			$ppStart=$start+$n; #bed format start is 0-based(included)
			$pp{"$chr:$ppStart-"}+=$score;
			$pp_seq{"$chr:$ppStart-"}=$len;
				#push @{$pp_seq{"$chr:$ppStart-"}},$len;
		
			}
			else 
			{ 
				$ppStart=$end-$n; #bed format end is 1-based(not included)
				$pp{"$chr:$ppStart+"}+=$score;
				$pp_seq{"$chr:$ppStart+"}=$len;
							#push @{$pp_seq{"$chr:$ppStart+"}},$len;
			}

		}
		close(IN);
	 
			
		open IN, $file2 or die "Cannot open $file2: $!\n";
		while(my $record=<IN>) 
		{
			chomp $record;
			my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
			my $len=0;
			$len=$end-$start;
			next if ($len>29 || $len<23);

			if ($strand eq '+') 
			{
					
				$pos{"$chr:$start+"}+=$score;
				$pos_seq{"$chr:$start+"}=$len;
			}
			else 
			{ 
				$pos{"$chr:$end-"}+=$score;
				$pos_seq{"$chr:$end-"}=$len;
			}	
		}
		close(IN);

	
     	#calculate the PPscore for each pair of piRNAs with 5'-5' distance n	
		foreach my $cor (keys %pos)
		{
			if ($cor && exists $pp{$cor})
			{
				$score{$n}+=$pos{$cor}*$pp{$cor} ; #product of reads count is the ppscore
				print PPSEQ "$cor\t$pos_seq{$cor}\t$pos{$cor}\t$pp_seq{$cor}\t$pp{$cor}\n" if ($n==10);
				#PPSEQ columns
				#chr:start(+/-) of the piRNA partner 
				#the length of piRNA partners (only record one length)
				#the abundance of piRNA partners
				#the length of piRNAs (only record one length)
				#the abundance of piRNA

			}
		}
		$score{$n}=0 if (!exists $score{$n});
		print PPSCORE "$n\t$score{$n}\n";
	
	} #n
	return \%score;	
