 #!/usr/bin/perl
 
use File::Basename;
 
my $inFile=$parameters->{input};
my $filename=fileparse($inFile);
my $outDir=$parameters->{outdir};
open IN, $inFile or die "could not find $inFile: $!";
open OUT, "> $outDir/$filename.baseFraction.txt"; 

my %A=();
my %C=(); 
my %G=(); 
my %T=();
my %count=();
my $lastPos=$parameters->{lenghth}-1; 
while(my $line=<IN>)
{
	chomp $line;
	my @l=split(/\t/,$line);
	my $piSeq=$l[$parameters->{seq}-1];
	my $readsNum=$l[$parameters->{readsNum}-1];
    
	
    foreach my $pos (0..$lastPos)
    {
        $A{$pos}{$piSeq}+=$readsNum if ( substr($piSeq,$pos,1) eq 'A');
        $C{$pos}{$piSeq}+=$readsNum if ( substr($piSeq,$pos,1) eq 'C');
        $G{$pos}{$piSeq}+=$readsNum if ( substr($piSeq,$pos,1) eq 'G');
     	$T{$pos}{$piSeq}+=$readsNum if ( substr($piSeq,$pos,1) eq 'T');
    }

}
my $countLastPos=$lastPos;

foreach my $pos (0..$countLastPos)
{
	my $aSpec=0;
	my $cSpec=0;
	my $gSpec=0;
	my $tSpec=0;
	my $totalSpec=0;
	
	$aSpec=scalar (keys %{$A{$pos}});
	$cSpec=scalar (keys %{$C{$pos}});
	$gSpec=scalar (keys %{$G{$pos}});
	$tSpec=scalar (keys %{$T{$pos}});
	
	$totalSpec=$aSpec+$cSpec+$gSpec+$gSpec;
	
		
	my $aFrac=0;
	my $cFrac=0;
	my $gFrac=0;
	my $tFrac=0;
	if($totalSpec != 0)
	{
		$aFrac=$aSpec/$totalSpec;
		$cFrac=$cSpec/$totalSpec;
		$gFrac=$gSpec/$totalSpec;
		$tFrac=$tSpec/$totalSpec;
	}

	print OUT "species\t$pos\t$aFrac\t$cFrac\t$gFrac\t$tFrac\t$totalSpec\n";
	
}


foreach my $pos (0..$countLastPos)
{
	
	my $aRead=0;
	my $cRead=0;
	my $gRead=0;
	my $tRead=0;
	my $totalRead=0;
	
	map {$aRead+=$_ } values %{$A{$pos}};
	map {$cRead+=$_ } values %{$C{$pos}};
	map {$gRead+=$_ } values %{$G{$pos}};
	map {$tRead+=$_ } values %{$T{$pos}};
	$totalRead=$aRead+$cRead+$gRead+$tRead;
	
	my $aFrac=0;
	my $cFrac=0;
	my $gFrac=0;
	my $tFrac=0;
	if($totalRead != 0)
	{
		$aFrac=$aRead/$totalRead;
		$cFrac=$cRead/$totalRead;
		$gFrac=$gRead/$totalRead;
		$tFrac=$tRead/$totalRead;
	}

	print OUT "reads\t$pos\t$aFrac\t$cFrac\t$gFrac\t$tFrac\t$totalRead\n";
}



  
  		sub usage
		{
		        print "\nUsage:$0\n\n\t";
		        print "REQUIRED\n\t";
		        print "-i  <inputfile>\n\t";
				print "-o  <outputdir>\n\t";				
				print "-p  <the column # of piRNA species>\n\t";
				print "-r  <the column # of reads number>\n\t";
				print "-l  <the maximal length to count>\n\t";
		        print "This perl script is count the frequency of nucleotide for each position of small RNAs\n";
				print "It's maintained by WEI WANG. If you have any questions, please contact wei.wang2\@umassmed.edu\n";
		        exit(1);
		}
		sub parse_command_line {
		        my($parameters, @ARGV) = @_;
		        my $next_arg;
		        while(scalar @ARGV > 0){
		                $next_arg = shift(@ARGV);
		                if($next_arg eq "-i"){ $parameters->{input} = shift(@ARGV); }		                
		                elsif($next_arg eq "-o"){ $parameters->{outdir} = shift(@ARGV); }		                
						elsif($next_arg eq "-p"){ $parameters->{seq}= shift(@ARGV); }
						elsif($next_arg eq "-r"){ $parameters->{readsNum} = shift(@ARGV); }
						elsif($next_arg eq "-l"){ $parameters->{lenghth}= shift(@ARGV); }
		                else{ print "Invalid argument: $next_arg"; usage(); }
		        }
		}		
  