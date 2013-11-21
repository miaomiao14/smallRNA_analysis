#!/usr/bin/perl
#08/01/2013
#WEI WANG
#supposed to be the publised version?

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;
use threads;
use Thread::Queue 3.02 qw( );


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




for (my $i=0; $i<@inputfiles; $i++) 
{
	for (my $j=0; $j<@inputfiles; $j++) 
	{
		
		my $file1="";	#$file[i] is the target strand,$file[j] is the guide strand
		$file1=fileparse($inputfiles[$i]);
		$file1 =~ /(.*)\.bed*.*/;
		my $filename1="";
		$filename1=$1;
    
		my $file2=""; #file2 is the guide strand
		$file2=fileparse($inputfiles[$j]);
		$file2 =~ /(.*)\.bed*.*/;
		my $filename2="";
		$filename2=$1;
   		
   		open PPZ, ">$outdir/${filename1}_${filename2}.zscore";
   		
		if($parameters->{trn})
		{
			my $Transposon=$parameters->{trn};
			print PPZ "Transposon\tPINGPONGFILE\tZSCORE\tPPSCORE10\n"; 
			print PPZ "$Transposon\t$filename1-$filename2\t";
			open PPSEQ, ">$outdir/${filename1}_$filename2.$Transposon.ppseq";
			open PPSCORE, ">$outdir/${filename1}_$filename2.$Transposon.ppscore";
		}
		else
		{
			print PPZ "PINGPONGFILE\tZSCORE\tPPSCORE10\n";
			print PPZ "$filename1-$filename2\t";
			open PPSEQ, ">$outdir/${filename1}_$filename2.ppseq";
			open PPSCORE, ">$outdir/${filename1}_$filename2.ppscore";
		}
 					#chr:start(+/-) of the piRNA partner 
					#the length of piRNA partners (only record one length)
					#the abundance of piRNA partners
					#the length of piRNAs (only record one length)
					#the abundance of piRNA   
    	print PPSEQ "POSOFPARTNER\tPARTNERLEN\tPARTNERABUNDANCE\tpiRNALEN\tpiRNAABUDNACE\n";
    	print PPSCORE "OVERLAP\tPPSCORE\n";
		my $X=0; my $Z=0; my %score=();
	
		foreach (my $n=1;$n<=20;$n++) #iterate with different 5'-5' distance of piRNAs from oposite strand
		{
			my %pp=(); my %pos=(); my %pp_seq=();my %pos_seq=();
			#read zip input file
			#the target strand
			if($file1=~/gz/)
			{
				my $gz="";
				$gz = gzopen($inputfiles[$i], "rb") or die "Cannot open $inputfiles[$i]: $gzerrno\n" ;
				while($gz->gzreadline(my $record) > 0)
				{ 
					chomp $record;
					my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
					my $len=0;
					$len=$end-$start;
					next if ($len>29 || $len<23); #the length range of piRNAs
					my $ppStart=0;
					if($parameters->{format} eq "bed")
					{
						if ($strand eq '+')  #if the piRNA is on + strand 
						{
							$ppStart=$start+$n; #bed format start is 0-based(included), get the start position of its piRNA partner
							$pp{"$chr:$ppStart-"}+=1; ##assuming each sequence(with n reads) will have n rows
													  ##record the chr:start position of its piRNA partner(-) as the key
							$pp_seq{"$chr:$ppStart-"}=$len;
							#push @{$pp_seq{"$chr:$ppStart-"}},$len; ##record the length of the piRNAs, not the length of its partner
						}
						else #if the piRNA is on - strand 
						{ 
							$ppStart=$end-$n; #bed format end is 1-based(not included)
							$pp{"$chr:$ppStart+"}+=1; ##record the chr:start position of its piRNA partner(+) as the key
							$pp_seq{"$chr:$ppStart+"}=$len;
							#push @{$pp_seq{"$chr:$ppStart+"}},$len;##record the length of the piRNAs, not the length of its partner
		
						}
					}
					if($parameters->{format} eq "bedscore")
					{
						if ($strand eq '+')  #if the piRNA is on + strand 
						{
							$ppStart=$start+$n; #bed format start is 0-based(included), get the start position of its piRNA partner
							$pp{"$chr:$ppStart-"}+=$score; ##assuming each sequence(with n reads) will have n rows
													  ##record the chr:start position of its piRNA partner(-) as the key
							$pp_seq{"$chr:$ppStart-"}=$len;
							#push @{$pp_seq{"$chr:$ppStart-"}},$len; ##record the length of the piRNAs, not the length of its partner
						}
						else #if the piRNA is on - strand 
						{ 
							$ppStart=$end-$n; #bed format end is 1-based(not included)
							$pp{"$chr:$ppStart+"}+=$score; ##record the chr:start position of its piRNA partner(+) as the key
							$pp_seq{"$chr:$ppStart+"}=$len;
							#push @{$pp_seq{"$chr:$ppStart+"}},$len;##record the length of the piRNAs, not the length of its partner
		
						}
					}
					
			     }
			     $gz->gzclose();
			}
			#not zip input file
			else
			{
				open IN, $inputfiles[$i] or die "Cannot open $inputfiles[$i]: $!\n";
				while(my $record=<IN>) 
				{
					chomp $record;
					my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
					my $len=0;
					$len=$end-$start;
					next if ($len>29 || $len<23);
					my $ppStart=0;
					if($parameters->{format} eq "bed")
					{
						if ($strand eq '+') 
						{
							$ppStart=$start+$n; #bed format start is 0-based(included)
							$pp{"$chr:$ppStart-"}+=1;
							$pp_seq{"$chr:$ppStart-"}=$len;
							#push @{$pp_seq{"$chr:$ppStart-"}},$len;
		
						}
						else 
						{ 
							$ppStart=$end-$n; #bed format end is 1-based(not included)
							$pp{"$chr:$ppStart+"}+=1;
							$pp_seq{"$chr:$ppStart+"}=$len;
							#push @{$pp_seq{"$chr:$ppStart+"}},$len;
						}
					}
					if($parameters->{format} eq "bedscore")
					{
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
					
				}
				close(IN);
	     	}
			#read zip input file
			#the guide strand
			if($file2=~/gz/)
			{
				my $gz="";
				$gz = gzopen($inputfiles[$j], "rb") or die "Cannot open $inputfiles[$j]: $gzerrno\n" ;
				while($gz->gzreadline(my $record) > 0)
				{ 
					chomp $record;
					my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
					my $len=0;
					$len=$end-$start;
					next if ($len>29 || $len<23);
					if($parameters->{format} eq "bed")
					{
						if ($strand eq '+') 
						{
							
							$pos{"$chr:$start+"}+=1; #record the chr:start position of the piRNAs
							$pos_seq{"$chr:$start+"}=$len;#record the length of the piRNAs
						}
						else 
						{ 
							
							$pos{"$chr:$end-"}+=1;
							$pos_seq{"$chr:$end-"}=$len;
		
						}
					}
					if($parameters->{format} eq "bedscore")
					{
						if ($strand eq '+') 
						{
							
							$pos{"$chr:$start+"}+=$score; #record the chr:start position of the piRNAs
							$pos_seq{"$chr:$start+"}=$len;#record the length of the piRNAs
						}
						else 
						{ 
							
							$pos{"$chr:$end-"}+=$score;
							$pos_seq{"$chr:$end-"}=$len;
		
						}
					}
			     }
			     $gz->gzclose();
			}
			#not zip inputfile
			else
			{
				open IN, $inputfiles[$j] or die "Cannot open $inputfiles[$j]: $!\n";
				while(my $record=<IN>) 
				{
					chomp $record;
					my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
					my $len=0;
					$len=$end-$start;
					next if ($len>29 || $len<23);
					if($parameters->{format} eq "bed")
					{
						if ($strand eq '+') 
						{
							
							$pos{"$chr:$start+"}+=1;
							$pos_seq{"$chr:$start+"}=$len;
						}
						else 
						{ 
							
							$pos{"$chr:$end-"}+=1;
							$pos_seq{"$chr:$end-"}=$len;
		
						}
					}
					if($parameters->{format} eq "bedscore")
					{
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
					
				}
				close(IN);
	     	}
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
