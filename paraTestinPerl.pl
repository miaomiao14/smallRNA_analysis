#!/usr/bin/perl
#11/29/2013
#WEI WANG
#supposed to be the publised version?

use strict;
use warnings;
use File::Basename;
use Compress::Zlib;

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
		`[ -s \$outidr/parafile.split.command ] && rm \$outidr/parafile.split.command;`;
		`echo -e "&split_chr($file1,$outdir)" >>\$outidr/parafile.split.command`;
		`echo -e "&split_chr($file2,$outdir)" >>\$outidr/parafile.split.command`;
		`ParaFly -c \$outidr/parafile.split.command -CPU 2 -failed_cmds \$outidr/parafile.split.failed_commands`;
   		
	}
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
