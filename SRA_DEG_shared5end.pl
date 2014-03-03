#!/usr/bin/perl
#03/03/2014
#WEI WANG
#check the shared 5 end between SRA and DEG

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;
BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "restrict_digts.pm";
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

		my $file1="";	#$file[i] is the target strand,$file[j] is the guide strand
		$file1=fileparse($inputfiles[0]);
		$file1 =~ /(.*)\.bed*.*/;
		my $filename1="";
		$filename1=$1;
    
		my $file2=""; #file2 is the guide strand
		$file2=fileparse($inputfiles[1]);
		$file2 =~ /(.*)\.bed*.*/;
		my $filename2="";
		$filename2=$1;
   		
   		open PPSCORE, ">$outdir/${filename1}_${filename2}.shared5end.stat.txt";
   		
		my %pp=(); my %pos=(); my %sharedpp=();my %sharedpos=();
			
			#read zip input file
			#the target strand
			if($file1=~/gz/)
			{
				my $gz="";
				$gz = gzopen($inputfiles[0], "rb") or die "Cannot open $inputfiles[0]: $gzerrno\n" ;
				while($gz->gzreadline(my $record) > 0)
				{ 
					chomp $record;
					my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
					my $len=0;
					$len=$end-$start;
					if($file1=~/SRA/)
					{
						next if ($len>29 || $len<23);
					}
					my $ppStart=0;
					if($parameters->{format} eq "bed")
					{							
							$pp{"$chr:$start,$strand"}+=1; ##assuming each sequence(with n reads) will have n rows
													  ##record the chr:start position of its piRNA partner(-) as the key
					}
					if($parameters->{format} eq "bedscore")
					{
							$pp{"$chr:$start,$strand"}+=$score; ##assuming each sequence(with n reads) will have n rows
					}
					if($parameters->{format} eq "normbed")
					{
							my ($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$record);
							$pp{"$chr:$bedstart,$strand"}+=$reads/$ntm;
					}
					
			     }
			     $gz->gzclose();
			}
			#not zip input file
			else
			{
				open IN, $inputfiles[0] or die "Cannot open $inputfiles[0]: $!\n";
				while(my $record=<IN>) 
				{
					chomp $record;
					my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
					my $len=0;
					$len=$end-$start;
					if($file1=~/SRA/)
					{
						next if ($len>29 || $len<23);
					}
					my $ppStart=0;
					if($parameters->{format} eq "bed")
					{

							$pp{"$chr:$start,$strand"}+=1;

					}
					if($parameters->{format} eq "bedscore")
					{

							$pp{"$chr:$start,$strand"}+=$score;

					}
					if($parameters->{format} eq "normbed")
					{
							my ($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$record);
							$pp{"$chr:$bedstart,$strand"}+=$reads/$ntm;
					}
					
				}
				close(IN);
	     	}
			#read zip input file
			#the guide strand
			if($file2=~/gz/)
			{
				my $gz="";
				$gz = gzopen($inputfiles[1], "rb") or die "Cannot open $inputfiles[1]: $gzerrno\n" ;
				while($gz->gzreadline(my $record) > 0)
				{ 
					chomp $record;
					my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
					my $len=0;
					$len=$end-$start;
					if($file2=~/SRA/)
					{
						next if ($len>29 || $len<23);
					}
					if($parameters->{format} eq "bed")
					{

							$pos{"$chr:$start,$strand"}+=1;

					}
					if($parameters->{format} eq "bedscore")
					{

							$pos{"$chr:$start,$strand"}+=$score;

					}
					if($parameters->{format} eq "normbed")
					{
							my ($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$record);
							$pos{"$chr:$bedstart,$strand"}+=$reads/$ntm;
					}
			     }
			     $gz->gzclose();
			}
			#not zip inputfile
			else
			{
				open IN, $inputfiles[1] or die "Cannot open $inputfiles[1]: $!\n";
				while(my $record=<IN>) 
				{
					chomp $record;
					my($chr,$start,$end,$name,$score,$strand) = split(/\t/,$record);
					my $len=0;
					$len=$end-$start;
					if($file2=~/SRA/)
					{
						next if ($len>29 || $len<23);
					}
					if($parameters->{format} eq "bed")
					{

							$pos{"$chr:$start,$strand"}+=1;

					}
					if($parameters->{format} eq "bedscore")
					{

							$pos{"$chr:$start,$strand"}+=$score;

					}
					if($parameters->{format} eq "normbed")
					{
							my ($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$record);
							$pos{"$chr:$bedstart,$strand"}+=$reads/$ntm;
					}
					
				}
				close(IN);
	     	}
	     	#calculate the PPscore for each pair of piRNAs with 5'-5' distance n	
			foreach my $cor (keys %pos)
			{
				if ($cor && exists $pp{$cor})
				{
					$sharedpp{$cor}+=$pp{$cor};
					$sharedpos{$cor}+=$pos{$cor};
				}
			}

			#pp
			my $totalppSpecies=0;
			my $totalppReads=0;
			$totalppSpecies=scalar (keys %pp);
			map { $totalppReads+=$_ } values %pp;
			
			my $sharedppSpecies=0;
			my $sharedppReads=0;
			$sharedppSpecies=scalar (keys %sharedpp);
			map { $sharedppReads+=$_ } values %sharedpp;
			
			$totalppReads=&restrict_num_decimal_digits($totalppReads,3);			
			$sharedppReads=&restrict_num_decimal_digits($sharedppReads,3);
			#fraction
			my $ppSpeciesF=0;
			$ppSpeciesF=$sharedppSpecies/$totalppSpecies;
			$ppSpeciesF=&restrict_num_decimal_digits($ppSpeciesF,3);
			my $ppReadsF=0;
			$ppReadsF=$sharedppReads/$totalppReads;
			$ppReadsF=&restrict_num_decimal_digits($ppReadsF,3);
			#pos
			my $totalposSpecies=0;
			my $totalposReads=0;
			$totalposSpecies=scalar (keys %pos);
			map { $totalposReads+=$_ } values %pos;
			
			my $sharedposSpecies=0;
			my $sharedposReads=0;
			$sharedposSpecies=scalar (keys %sharedpos);
			map { $sharedposReads+=$_ } values %sharedpos;
			#fraction
			my $posSpeciesF=0;
			$posSpeciesF=$sharedposSpecies/$totalposSpecies;
			$posSpeciesF=&restrict_num_decimal_digits($posSpeciesF,3);
			my $posReadsF=0;
			$posReadsF=$sharedposReads/$totalposReads;									
			$posReadsF=&restrict_num_decimal_digits($posReadsF,3);
			$totalposReads=&restrict_num_decimal_digits($totalposReads,3);			
			$sharedposReads=&restrict_num_decimal_digits($sharedposReads,3);
			
			
			
			print PPSCORE "$sharedppSpecies\t$totalppSpecies\t$ppSpeciesF\t$sharedppReads\t$totalppReads\t$ppReadsF\t$sharedposSpecies\t$totalposSpecies\t$posSpeciesF\t$sharedposReads\t$totalposReads\t$posReadsF\n";
			


close(PPSCORE);


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
	print "-f  <input format [bed|bedscore|normbed]>\n\t";
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
