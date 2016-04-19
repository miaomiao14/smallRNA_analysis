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

#pay spcial attention to hash2 because the sequences are not stored in degradome files

#01/22/2014
#wei.wang2@umassmed.edu

#05/15-2014
#change windowsize from 20 to 16

#08/07/2014
#modify to fit ghpcc

#09/18
#windowsize -50 to 50
#allow mismatch only at g1
#$fa="$HOMEDIR/smallRNApipeline/pipeline_dm/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
#$fa="/home/wangw1/data/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
$fa="/home/wangw1/data/common/dmel-all-chromosome-r5.5_TAS.fasta";
open IN, $fa;
while(my $line=<IN>) { if ($line=~/>(.+) type/) { $chr="chr$1";} else { chomp $line; $genome{$chr}=$line;}}

#/home/hanb/nearline/small_RNA_Pipeline/common_files/dm3/genome.fa

%hash1=();
%hash2=();
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
		#$hash2{$file}{substr($_[4],0,16)}+=$_[5]/$_[6]; #here is the problem for DEG
		
		if($l[3] eq '+')
		{
			$seqstart=$l[1]-1; #norm.bed is 1-based
            $hash2{$file}{substr($genome{$l[0]},$seqstart,16)}+=$l[5]/$l[6]; #hash2 is the guide
		}
		else
		{
			$seqstart=$l[2]-16;
     		$seqtemp=substr($genome{$l[0]},$seqstart,16);
     		$seqtemp=&revfa($seqtemp);
     		$hash2{$file}{$seqtemp}+=$l[5]/$l[6]; #hash2 is the guide
		}
		
		
		for ($n=-50;$n <=50;$n++)
		{
			if ($l[3] eq '+')
			{
				$start=$l[1]+$n-17;
				$str=substr($genome{$l[0]},$start,16);
				$str=&revfa($str);
                $hash1{$file}{$n}{$str}+=$l[5]/$l[6];     #str is the target; hash1
                push @{$hash3{$file}{$n}{$str}},$l[4] if($n==10); #hash3, key is the target, value is the guide
                
                
                
			}
			else 
			{
				$start=$l[2]-$n;
				$str=substr($genome{$l[0]},$start,16);
      			$hash1{$file}{$n}{$str}+=$l[5]/$l[6]; #str is the target; hash1
     			push @{$hash3{$file}{$n}{$str}},$l[4] if($n==10); #hash3, key is the target, value is the guide
     			
     			

      		}
		}
	}
	if ($total{$file}>10) 
	{
		$seqFile="$OUTDIR/$file.seq";
		open OUT, ">$seqFile";
		foreach my $prefix (keys %{$hash2{$file}}) { print OUT "$prefix\t$hash2{$file}{$prefix}\n" if (length($prefix)==16)}#guide as query
		close(OUT);
		for ($n=1;$n <=16;$n++) 
		{
		$fa="$OUTDIR/$file.ref.$n.fa";
	    $indexb="$OUTDIR/$file.$n";
		open OUT, ">$fa";
		foreach my $prefix (keys %{$hash1{$file}{$n}}) { print OUT ">$prefix\t$hash1{$file}{$n}{$prefix}\n$prefix\n" if (length($prefix)==16)} #target as index
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
					$score{$n}+=$hash1{$file1}{$n}{$l[1]}*$hash2{$file2}{$l[2]}/$NTM{$l[2]};
					if ($n==10)
		     		{
		      			print PPSEQ "$l[2]\t$hash2{$file2}{$l[2]}\t$hash1{$file1}{$n}{$l[1]}\t" ; #hash2 is the guide; 
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


#`rm $OUTDIR/*.ebwt`;
#`rm $OUTDIR/*.bowtie.out`;
#`rm $OUTDIR/*.seq`;