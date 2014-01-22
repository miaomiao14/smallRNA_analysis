#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
use Compress::Zlib;
# rule is p1 and p17-21 doesn't need to pair but p2-10 need, the rest 1mm
# simplifized version with prefix 16nt
# input as norm.bed ( no header)


open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while(<IN>) { if (/>(.+) type/) { $chr="chr$1";} else { chomp; $genome{$chr}=$_;}}

%hash1=();
%hash2=();
@hash3=();
$OUTDIR=$ARGV[3];
# make bowtie index
for ($i=0; $i<$ARGV[2]; $i++)
{
	my $file=fileparse($ARGV[$i]);
	@namefield=split(/\./,$file);
	$name=$namefield[2]."_".$namefield[1]."_".$namefield[7];

	$file=$name;
	open IN, $ARGV[$i];
	while(<IN>) 
	{
		chomp; 
		split(/\t/);  
		#next if (length($_[4])>29 || length($_[4])<23); #can not use it for DEG
		next if (/data/);
		$total{$file}+=$_[5]/$_[6];
		#$hash2{$file}{substr($_[4],0,16)}+=$_[5]/$_[6]; #here is the problem for DEG
		
		if($_[3] eq '+')
		{
			$seqstart=$_[1]-1; #norm.bed is 1-based
            $hash2{$file}{substr($genome{$_[0]},$seqstart,16)}+=$_[5]/$_[6]; #hash2 is the guide
		}
		else
		{
			$seqstart=$_[2]-16;
     		$seqtemp=substr($genome{$_[0]},$seqstart,16);
     		$seqtemp=&revfa($seqtemp);
     		$hash2{$file}{$seqtemp}+=$_[5]/$_[6]; #hash2 is the guide
		}
		
		
		for ($n=1;$n<=20;$n++)
		{
			if ($_[3] eq '+')
			{
				$start=$_[1]+$n-17;
				$str=substr($genome{$_[0]},$start,16);
				$str=&revfa($str);
                $hash1{$file}{$n}{$str}+=$_[5]/$_[6];     #str is the target; hash1
                push @{$hash3{$file}{$n}{$str}},$_[4] if($n==10); #hash3, key is the target, value is the guide
                
                
                
			}
			else 
			{
				$start=$_[2]-$n;
				$str=substr($genome{$_[0]},$start,16);
      			$hash1{$file}{$n}{$str}+=$_[5]/$_[6]; #str is the target; hash1
     			push @{$hash3{$file}{$n}{$str}},$_[4] if($n==10); #hash3, key is the target, value is the guide
     			
     			

      		}
		}
	}
	if ($total{$file}>10) 
	{
		$seqFile="$OUTDIR/$file.seq";
		open OUT, ">$seqFile";
		foreach (keys %{$hash2{$file}}) { print OUT "$_\t$hash2{$file}{$_}\n" if (length($_)==16)}#guide as query
		close(OUT);
		for ($n=1;$n<=20;$n++) 
		{
		$fa="$OUTDIR/$file.ref.$n.fa";
	    $indexb="$OUTDIR/$file.$n";
		open OUT, ">$fa";
		foreach (keys %{$hash1{$file}{$n}}) { print OUT ">$_\t$hash1{$file}{$n}{$_}\n$_\n" if (length($_)==16)} #target as index
		`bowtie-build $fa $indexb && rm $fa`;

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
			foreach ($n=1;$n<=20;$n++) 
			{
				# file1 as ref; target as ref;
		 		`bowtie $OUTDIR/$file1.$n -r -a -v 1 -p 8 $OUTDIR/$file2.seq --suppress 1,4,6,7 | grep + > $OUTDIR/$file2.$file1.$n.bowtie.out`;
		  		%NTM=();
		  		open IN, "$OUTDIR/$file2.$file1.$n.bowtie.out";
		  		
				while(<IN>)
				{
					chomp; split(/\t/);  
					if ($_[3]=~/(\d+):/) { next if ($1<=9 && $1>=1);}  # no seed mismatches 0 based
					$NTM{$_[2]}++; #a guide can map to multiple targets
		  		}
		  		close(IN);
		  		open IN, "$OUTDIR/$file2.$file1.$n.bowtie.out";
				while(<IN>)
				{
					chomp; split(/\t/);
					if ($_[3]=~/(\d+):/){next if ($1<=9 && $1>=1);}  # no seed mismatches 0 based
					$score{$n}+=$hash1{$file1}{$n}{$_[1]}*$hash2{$file2}{$_[2]}/$NTM{$_[2]};
					if ($n==10)
		     		{
		      			print PPSEQ "$_[2]\t$hash2{$file2}{$_[2]}\t$hash1{$file1}{$n}{$_[1]}\t" ; #hash2 is the guide; 
		      			local $"=',';
		      			print PPSEQ "@{$hash3{$file1}{$n}{$_[1]}}";#hash3, key is the target seq, value is the guide seq
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
	 	
	 	$N1=`match.pl $OUTDIR/$file2.$file1.ppseq $OUTDIR/$file2.seq | sumcol+ 2`; chomp($N1);
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


`rm $OUTDIR/*.ebwt`;
`rm $OUTDIR/*.bowtie.out`;
`rm $OUTDIR/*.seq`;