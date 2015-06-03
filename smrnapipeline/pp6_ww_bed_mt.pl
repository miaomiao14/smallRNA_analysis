#!/usr/bin/perl
#11/29/2013
#WEI WANG
#supposed to be the publised version?

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;
$ENV{PATH} = "/home/wangw1/git/smallRNA_analysis/Utils:$ENV{PATH}";
my $parameters={};

if(scalar(@ARGV)==0)
{
	usage();
	}

&parse_command_line($parameters, @ARGV);

my $inputfile1=fileparse($parameters->{in_file1});
my @inputfiles=();
push @inputfiles, $parameters->{in_file1};
my $inputfile2="";

my $outdir=$parameters->{outdir};
my $format=$parameters->{format};

if($parameters->{in_file2})
{
	$inputfile2=fileparse($parameters->{in_file2});
	push @inputfiles, $parameters->{in_file2};
	
}



#main
for (my $i=0; $i<@inputfiles; $i++) 
{
	for (my $j=0; $j<@inputfiles; $j++) 
	{
		
		my $file1="";
		$file1=$inputfiles[$i];	#$file[i] is the target strand,$file[j] is the guide strand
		my $filetemp1="";
		$filetemp1=fileparse($inputfiles[$i]);
		$filetemp1 =~ /(.*)\.bed*.*/;
		my $filename1="";
		$filename1=$1;
    
		my $file2="";
		$file2=$inputfiles[$j]; #file2 is the guide strand
		my $filetemp2="";
		$filetemp2=fileparse($inputfiles[$j]);
		$filetemp2 =~ /(.*)\.bed*.*/;
		my $filename2="";
		$filename2=$1;
		my $outdir=$parameters->{outdir};
		`[ -s $parameters->{outdir}/parafile.split.command ] && rm $parameters->{outdir}/parafile.split.command;`;
		`echo -e " split_chr.pl $file1 $parameters->{outdir}" >> $parameters->{outdir}/parafile.split.command`;
		if($filename2 ne $filename1)
		{
		`echo -e " split_chr.pl $file2 $parameters->{outdir}" >> $parameters->{outdir}/parafile.split.command`;
		}
		
		`ParaFly -c $parameters->{outdir}/parafile.split.command -CPU 2 -failed_cmds $parameters->{outdir}/parafile.split.failed_commands`;
		my %chrs1=();my %chrs2=();
   		open CHR,"$parameters->{outdir}/$filetemp1.chrlist.txt" or die "Cannot open $parameters->{outdir}/$filetemp1.chrlist.txt: $!\n"; #the output from split_chr.pl
   		while(my $line=<CHR>)
   		{
   			chomp $line;
   			my @l=split(/\t/,$line);
   			$chrs1{$l[0]}=1;
   		}
   		close(CHR);
   		
   		open CHR,"$parameters->{outdir}/$filetemp2.chrlist.txt" or die "Cannot open $parameters->{outdir}/$filetemp2.chrlist.txt: $!\n"; #the output from split_chr.pl
   		while(my $line=<CHR>)
   		{
   			chomp $line;
   			my @l=split(/\t/,$line);
   			$chrs2{$l[0]}=1;
   		}
   		close(CHR);
   		
   		open PPZ, ">$outdir/${filename1}.${filename2}.zscore";
   		
		if($parameters->{trn})
		{
			my $Transposon=$parameters->{trn};
			print PPZ "Transposon\tPINGPONGFILE\tZSCORE\tPPSCORE10\n"; 
			print PPZ "$Transposon\t$filename1-$filename2\t";
			open PPSEQ, ">$outdir/${filename1}.$filename2.$Transposon.ppseq";
			open PPSCORE, ">$outdir/${filename1}.$filename2.$Transposon.ppscore";
		}
		else
		{
			print PPZ "PINGPONGFILE\tZSCORE\tPPSCORE10\n";
			print PPZ "$filename1-$filename2\t";
			open PPSEQ, ">$outdir/${filename1}.$filename2.ppseq";
			open PPSCORE, ">$outdir/${filename1}.$filename2.ppscore";
		}
 		#chr:start(+/-) of the piRNA partner 
		#the length of piRNA partners (only record one length)
		#the abundance of piRNA partners
		#the length of piRNAs (only record one length)
		#the abundance of piRNA   
    	print PPSEQ "POSOFPARTNER\tPARTNERLEN\tPARTNERABUNDANCE\tpiRNALEN\tpiRNAABUDNACE\n";
    	print PPSCORE "OVERLAP\tPPSCORE\n";
		my $X=0; my $Z=0;
		my %score=(); 
		my $scoreref="";
		my %scorechr=();
		if($parameters->{format} eq "bed")
		{
			foreach my $chromosome (%chrs1)
			{
				if ($chrs2{$chromosome})
				{
					`[ -s $parameters->{outdir}/parafile.ppscore.command ] && rm $parameters->{outdir}/parafile.ppscore.command;`;
					#`echo -e " split_chr.pl $file1 $parameters->{outdir}" >> $parameters->{outdir}/parafile.split.command`;
					`echo -e " ppscore_calculation.pl $outdir/$filetemp1.$chromosome $outdir/$filetemp2.$chromosome bed $chromosome" >> $parameters->{outdir}/parafile.ppscore.command `;
				}
			}

		} #format bed
		if($parameters->{format} eq "bedscore")
		{
			foreach my $chromosome (%chrs1)
			{
				if ($chrs2{$chromosome})
				{
					`[ -s $parameters->{outdir}/parafile.ppscore.command ] && rm $parameters->{outdir}/parafile.ppscore.command;`;
					#`echo -e " split_chr.pl $file1 $parameters->{outdir}" >> $parameters->{outdir}/parafile.split.command`;
					`echo -e " ppscore_calculation.pl $outdir/$filetemp1.$chromosome $outdir/$filetemp2.$chromosome bed $chromosome" >> $parameters->{outdir}/parafile.ppscore.command `;
				}
			}
		} #format bedscore
		
		my $CPUN=`wc -l $parameters->{outdir}/parafile.ppscore.commands|cut -f1 -d" "`;
		$CPUN=int($CPUN/4);
		`ParaFly -c $parameters->{outdir}/parafile.ppscore.command -CPU $CPUN -failed_cmds $parameters->{outdir}/parafile.split.failed_commands`;
		
		#sum the score from each chromosomes
		for(my $n=1; $n<21; $n++)
		{
			foreach my $chromosome (%chrs1)
			{
				if ($chrs2{$chromosome})
				{
					open PPSCORE "$outdir/";
					$score{$n}+=$scorechr{$chromosome}->{$n};
				}
			}
		}
		
		
		if ($score{10}) { $X=$score{10}; delete $score{10};} #mean and standard deviation exclude ppscore10
		my $std=0;
		$std=&standard_deviation(values %score);
		if ((scalar %score)>9 && $std>0) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
		printf PPZ '%.2f', $Z;
		print PPZ "\t$X\n";
	}

}
close(PPSCORE);
close(PPSEQ);
close(PPZ);

sub mean 
{
	my $count=0;
	my(@numbers) =@_;
	foreach (@_) { $count+=$_;}
	return $count/(scalar @_);
}

sub standard_deviation 
{
	my(@numbers) = @_;
	#Prevent division by 0 error in case you get junk data
	return undef unless(scalar(@numbers));
	
	# Step 1, find the mean of the numbers
	my $total1 = 0;
	foreach my $num (@numbers) {
	$total1 += $num;
	}
	my $mean1 = $total1 / (scalar @numbers);
	
	# Step 2, find the mean of the squares of the differences
	# between each number and the mean
	my $total2 = 0;
	foreach my $num (@numbers) {
	$total2 += ($mean1-$num)**2;
	}
	my $mean2 = $total2 / (scalar @numbers);
	
	# Step 3, standard deviation is the square root of the
	# above mean
	my $std_dev = sqrt($mean2);
	return $std_dev;
}

sub usage
{
	print "\nUsage:$0\n\n\t";
	print "REQUIRED\n\t";
	print "-i  <inputfile1>\n\t";
	print "-j  <inputfile2 [optional]>\n\t";
	print "-o  <outputdir>\n\t";
	print "-f  <input format [bed|bedscore]>\n\t";
	print "-t  <transposon name[optional]>\n\t";
	print "This perl script is to calculate the Ping-Pong Zscore of smallRNAs within length range 23-29nt\n\t";
	print "The input is bed file which includes both +(or sense) and -(or antisense) mappers\n\t";
	print "It's maintained by WEI WANG. If you have any questions, please contact wei.wang2\@umassmed.edu\n";
	
	exit(1);
	}
sub parse_command_line {
        my($parameters, @ARGV) = @_;
        my $next_arg;
        while(scalar @ARGV > 0){
                $next_arg = shift(@ARGV);
                if($next_arg eq "-i"){ $parameters->{in_file1} = shift(@ARGV); }
                elsif($next_arg eq "-j"){ $parameters->{in_file2} = shift(@ARGV); }
                elsif($next_arg eq "-o"){ $parameters->{outdir} = shift(@ARGV); }
                elsif($next_arg eq "-f"){ $parameters->{format} = shift(@ARGV); }
                elsif($next_arg eq "-t"){ $parameters->{trn} = shift(@ARGV); }
                else{ print "Invalid argument: $next_arg"; usage(); }
        }
}

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


