#!/usr/bin/perl


BEGIN { unshift @INC,"$ENV{\"HOME\"}/git/smallRNA_analysis/Utils/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
use Compress::Zlib;
# rule is p1 and p17-21 doesn't need to pair but p2-10 need, the rest 1mm
# simplifized version with prefix 16nt
# input as norm.bed ( no header)

#adjust the script for PP8 between smRNAs and Degradome
#input are norm.bed file 

#pay spcial attention to readMap2 because the sequences are not stored in degradome files

#01/22/2014
#wei.wang2@umassmed.edu

#05/15-2014
#change windowsize from 20 to 16

#08/07/2014
#modify to fit ghpcc

#09/18/2014
#windowsize -50 to 50
#allow mismatch only at g1
#$fa="$HOMEDIR/smallRNApipeline/pipeline_dm/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
#$fa="/home/wangw1/data/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";

#09/24/2014
#format

my $numArgs = $#ARGV + 1;
my $format='normbed';
my $CLONE='SRA';
my $winSize=100;
my $inFileName1='';
my $inFileName2='';
my $end1='';
my $end2='';
my $outDir=getcwd;
my $help="Help";
my $ok = GetOptions( 'help|?' => \$help,'inputFile1|i=s' => \$inFileName1, 'end1|e1=s' => \$end1, 'inputFile2|j=s' => \$inFileName2,'end2|e2=s' => \$end2, 'outputDir|o=s' => \$outDir, 'format|F=s' => \$format, 'Clone|C=s' => $CLONE, 'windowSize|w=i' => \$winSize );

if($numArgs < 5 || !$ok || !$inFileName1 || !$outDir )
{
  Usage();
  exit;
}


$fa="/home/wangw1/data/common/dmel-all-chromosome-r5.5_TAS.fasta";
open IN, $fa;
while(my $line=<IN>) { if ($line=~/>(.+) type/) { $chr="chr$1";} else { chomp $line; $genome{$chr}=$line;}}

#/home/hanb/nearline/small_RNA_Pipeline/common_files/dm3/genome.fa

my %readMap1=();
my %readMap2=();
@hash3=();
$OUTDIR=$ARGV[3];
# make bowtie index
for ($i=0; $i<$ARGV[2]; $i++)
{
	my $file=fileparse($ARGV[$i]);
	$filename=$file;
	@namefield=split(/\./,$file);
	$name=$namefield[2]."_".$namefield[1]."_".$namefield[7];

	$file=$name;
	
	parseInputFile(\%readMap1,$inFileName1,$end1);
	
	$gz = gzopen($ARGV[$i], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
	while($gz->gzreadline($line) > 0)
	#open IN, $ARGV[$i];
	#while(<IN>) 
	{
		chomp $line; 
		@l=split(/\t/,$line);
		if($filename=~/SRA/ or $filename=~/LRA/)
		{  
			next if (length($l[4])>29 || length($l[4])<23); #can not use it for DEG
		}
		next if ($line=~/data/);
		$total{$file}+=$l[5]/$l[6];
		#$readMap2{$file}{substr($_[4],0,16)}+=$_[5]/$_[6]; #here is the problem for DEG
		
		if($l[3] eq '+')
		{
			$seqstart=$l[1]-1; #norm.bed is 1-based
            $readMap2{$file}{substr($genome{$l[0]},$seqstart,16)}+=$l[5]/$l[6]; #readMap2 is the guide
		}
		else
		{
			$seqstart=$l[2]-16;
     		$seqtemp=substr($genome{$l[0]},$seqstart,16);
     		$seqtemp=&revfa($seqtemp);
     		$readMap2{$file}{$seqtemp}+=$l[5]/$l[6]; #readMap2 is the guide
		}
		
		
		for ($n=-50;$n <=50;$n++)
		{
			if ($l[3] eq '+')
			{
				$start=$l[1]+$n-17;
				$str=substr($genome{$l[0]},$start,16);
				$str=&revfa($str);
                $readMap1{$file}{$n}{$str}+=$l[5]/$l[6];     #str is the target; readMap1
                push @{$hash3{$file}{$n}{$str}},$l[4] if($n==10); #hash3, key is the target, value is the guide
                
                
                
			}
			else 
			{
				$start=$l[2]-$n;
				$str=substr($genome{$l[0]},$start,16);
      			$readMap1{$file}{$n}{$str}+=$l[5]/$l[6]; #str is the target; readMap1
     			push @{$hash3{$file}{$n}{$str}},$l[4] if($n==10); #hash3, key is the target, value is the guide
     			
     			

      		}
		}
	}
	if ($total{$file}>10) 
	{
		$seqFile="$OUTDIR/$file.seq";
		open OUT, ">$seqFile";
		foreach my $prefix (keys %{$readMap2{$file}}) { print OUT "$prefix\t$readMap2{$file}{$prefix}\n" if (length($prefix)==16)}#guide as query
		close(OUT);
		for ($n=1;$n <=16;$n++) 
		{
		$fa="$OUTDIR/$file.ref.$n.fa";
	    $indexb="$OUTDIR/$file.$n";
		open OUT, ">$fa";
		foreach my $prefix (keys %{$readMap1{$file}{$n}}) { print OUT ">$prefix\t$readMap1{$file}{$n}{$prefix}\n$prefix\n" if (length($prefix)==16)} #target as index
		`bowtie-build $fa $indexb && rm $fa`;
		#`bowtie-build $fa $indexb `;

	 	}
		close (OUT);
	}
}



open PPZ, ">$OUTDIR/zscore.out";

# bowtie mapping and score calculating
for ($i=0; $i<$ARGV[2]; $i++)
{
	for ($j=0; $j<$ARGV[2]; $j++)
	{
      
	$file1=fileparse($ARGV[$i]); 
	@namefield=split(/\./,$file1);
	$name1=$namefield[2]."_".$namefield[1]."_".$namefield[7];
	$file1=$name1;
	
	$file2=fileparse($ARGV[$j]);
	@namefield=split(/\./,$file2);
	$name2=$namefield[2]."_".$namefield[1]."_".$namefield[7];
	$file2=$name2;

	if ($total{$file1}<10 || $total{$file2}<10) { print "$file2-$file1\t-10\t"; print PPZ "$file2-$file1\t-10\t";} 
	else 
		{
	 		$X=0; $Z=0; %score=();$count_N=0;%scorenf=();
			
			open PPSCORE, ">$OUTDIR/$file2.$file1.pp";
			open PPSEQ, ">$OUTDIR/$file2.$file1.ppseq";
			foreach ($n=-50;$n <=50;$n++) 
			{
				# file1 as ref; target as ref;
		 		`bowtie $OUTDIR/$file1.$n -r -a -v 1 -p 8 $OUTDIR/$file2.seq --suppress 1,4,6,7 | grep + > $OUTDIR/$file2.$file1.$n.bowtie.out`;
		  		%NTM=();
		  		open IN, "$OUTDIR/$file2.$file1.$n.bowtie.out";
		  		
				while(my $line=<IN>)
				{
					chomp $line; @l= split(/\t/,$line);  
					if ($l[3]=~/(\d+):/) { next if ($1<=16 && $1>=1);}  # no seed mismatches 0 based
					$NTM{$l[2]}++; #a guide can map to multiple targets
		  		}
		  		close(IN);
		  		open IN, "$OUTDIR/$file2.$file1.$n.bowtie.out";
				while(my $line=<IN>)
				{
					chomp $line; @l=split(/\t/,$line);
					if ($l[3]=~/(\d+):/){next if ($1<=16 && $1>=1);}  # no seed mismatches 0 based
					$score{$n}+=$readMap1{$file1}{$n}{$l[1]}*$readMap2{$file2}{$l[2]}/$NTM{$l[2]};
					if ($n==10)
		     		{
		      			print PPSEQ "$l[2]\t$readMap2{$file2}{$l[2]}\t$readMap1{$file1}{$n}{$l[1]}\t" ; #readMap2 is the guide; 
		      			local $"=',';
		      			print PPSEQ "@{$hash3{$file1}{$n}{$l[1]}}";#hash3, key is the target seq, value is the guide seq
		      			print PPSEQ "\n";
		     		}
		  		}
		  		close(IN);
		  		  
		
		     	$score{$n}=0 if (!exists $score{$n});
		     	$scorenf{$n}=$score{$n}*1000000000000/$total{$file1}/$total{$file2};
		     	print PPSCORE "$n\t$score{$n}\t$scorenf{$n}\n";
		     	$count_N++ if ($score{$n}>0);
			}
			close(PPSEQ);
			close(PPSCORE);
			
			
		$X=$score{10}; delete $score{10};
		$std=&std(values %score); 
		if ($std>0 && $count_N>=5) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
		print PPZ "$file2-$file1\t$Z\t";
		print "$file2-$file1\t$Z\t";
	 	
	 	$N1=`/home/wangw1/bin/match.pl $OUTDIR/$file2.$file1.ppseq $OUTDIR/$file2.seq | sumcol+ 2`; chomp($N1);
	 	
		if ($Z!=-10) 
		{
			print PPZ "$X\t$N1\t$total{$file2}\t",$N1/$total{$file2},"\t","$total{$file1}\t",$X*1000000000000/$total{$file1}/$total{$file2},"\n";
			print "$X\t$N1\t$total{$file2}\t",$N1/$total{$file2},"\t","$total{$file1}\t",$X*1000000000000/$total{$file1}/$total{$file2},"\n";
		} 
		else
		{
			print PPZ "$X\tNA\t$total{$file2}\tNA\t$total{$file1}\tNA\n";
			print "$X\tNA\t$total{$file2}\tNA\t$total{$file1}\tNA\n";
		}
		
		}#else

	}#j

}#i
close(PPZ);


sub parseInputFile
{
  		my ($readMapRef, $file, $anchor) = @_;
		$anchor=5;
		my $fileHandle = FileHandle->new;
		if ($fileHandle->open("< $file"))
		{

        	if($format eq "normbed")
        	{
        		while(my $line = $fileHandle->getline())
				{ 
					next if ($line=~/data/ or $line=~/\#/);
					chomp $line; 
					$line=~s/\s+/\t/g; #change all spaces to tab
					
					my $chr;
					my $bedstart;
					my $bedend;
					my $strand;
					my $seq;
					my $reads;
					my $ntm;
					my $len;
        		
	        		($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$line);
					$bedstart=$bedstart-1;
					
					
					&tagRecord($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm,$readMapRef,$anchor);
					
					
				}#while
        	}
			if($format eq "bed")
			{
				while(my $line = $fileHandle->getline())
				{ 
					next if ($line=~/data/ or $line=~/\#/);
					chomp $line; 
					$line=~s/\s+/\t/g; #change all spaces to tab
					
					my $chr;
					my $bedstart;
					my $bedend;
					my $strand;
					my $seq;
					my $reads;
					my $ntm;
					my $len;
					my $ntmreads;
					
					($chr,$bedstart,$bedend,$seq,$ntmreads,$strand)= split(/\t/,$line);
					$reads=$ntmreads;
					$ntm=1;
					
					
					&tagRecord($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm,$readMapRef,$anchor);
				}#while
        	}
        	if($format eq "bed2")#bo's definition
			{
				
				while(my $line = $fileHandle->getline())
				{ 
					next if ($line=~/data/ or $line=~/\#/);
					chomp $line; 
					$line=~s/\s+/\t/g; #change all spaces to tab
					
					my $chr;
					my $bedstart;
					my $bedend;
					my $strand;
					my $seq;
					my $reads;
					my $ntm;
					my $len;
					($chr,$bedstart,$bedend,$reads,$ntm,$strand,$seq)= split(/\t/,$line);
					
	
					&tagRecord($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm,$readMapRef,$anchor);
        	}#while		
        			
		}#if
	}#filehandle
	$fileHandle->close();
} #sub
		

sub tagRecord					 
{
	my ($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm,$readMapRef,$anchor) = @_ ;
	my $len=$bedend-$bedstart;
	if($CLONE eq "SRA") #length range to filter for piRNAs; for DEG, no filtering should be processed
	{
		return if ($len>29 or $len<23);
	}
	if($anchor==5)
	{
		my $nuc=substr($seq,0,1);
		
		if ($strand eq '+')
		{
        	$readMapRef->{$strand}{$chr}{$bedstart}->[1]+=$reads/$ntm; 
			$readMapRef->{$strand}{$chr}{$bedstart}->[0]=$nuc;
        	#$plus_end{$chr}{$bedend}=$l[2];
    	}
    	else
    	{
			my $fiveEnd=$bedend-1;
        	$readMapRef->{$strand}{$chr}{$fiveEnd}->[1]+=$reads/$ntm;
			$readMapRef->{$strand}{$chr}{$fiveEnd}->[0]=$nuc;
    	}
	}
	if($anchor==3)
	{
		my $nuc=substr($seq,-1);
		if ($strand eq '+')
		{
			#this 3' end is closed
			my $threeEnd=$bedend-1;
        	$readMapRef->{$strand}{$chr}{$threeEnd}->[1]+=$reads/$ntm; 
			$readMapRef->{$strand}{$chr}{$threeEnd}->[0]=$nuc;
        	#$plus_end{$chr}{$bedend}=$l[2];
    	}
    	else
    	{
			#if minus strand, 
			my $threeEnd=$bedstart;
        	$readMapRef->{$strand}{$chr}{$threeEnd}->[1]+=$reads/$ntm;
			$readMapRef->{$strand}{$chr}{$threeEnd}->[0]=$nuc;
    	}
	}
}#sub




sub mean 
{
	my $count=0;
	my(@numbers) =@_;
	foreach my $num (@numbers) { $count+=$num;}
	return $count/(scalar @numbers);
}

sub std 
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

sub Usage
{
	print "\n=======================================================================\n"; 
	print "\nUSAGE:$0\n\n -i <inputFile1> -j <inputFile2> -e1 <file1 End of Sequences> -e2 <file2 End of Sequences> -o <outputDir> -F <bed|normbed|bed2> -C <SRA|DEG> -w <windowSize of distance>\n\n"; 
	print " [-i <file>]\tinput file1 containing tags/reads in normbed/BED/bed2 format\n";
	print " [-j <file>]\tinput file2 containing tags/reads in normbed/BED/bed2 format\n";
	print " [-e1 <str>]\twhich end[5|3] of input file1\n";
	print " [-e2 <str>]\twhich end[5|3] of input file2\n";
	print " [-o <file>]\toutput file into which directory will be stored (default: current working directory)\n";
	print " [-F <str>]\tinput-file type, currently bed,bed2(customized) and normbed are supported\n\t\t(default: normbed)\n";
	print " [-C <str>]\tcloning method (SRA or DEG, length matters, default:SRA)\n";
	print " [-w <int>]\tscanning window size downstread of each signal (default: 100)\n";
 	print "This perl script is count the frequency of 5'-5'end distances of smallRNAs(23-29) from the same strand\n";
	print  "\n=======================================================================\n\n";
	exit(1); 
} 


#`rm $OUTDIR/*.ebwt`;
#`rm $OUTDIR/*.bowtie.out`;
#`rm $OUTDIR/*.seq`;