#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
use Compress::Zlib;
# rule is p1 and p17-21 doesn't need to pair but p2-16 need
# simplifized version with prefix 16nt
# input as norm.bed ( no header)


##This is another version of pp8,only allowing for 1mm at the 1st position
##calculating the frequency of VA,HC,DG,BT
##V is not equal to U
##V represents the first sequence composition of guide strand
##A represents the 10th sequence composition of target strand
## HC,DG,BT as controls
## separate Ago3-Aub and Aub-Ago3 PP

#11/05/2013
#fasta fiel as a input argument
#outdir added
#indexflag to indicate if we need to rebuild the index or not

#02/10/2014
#Phil wants to separate the cis-1U10A from trans-1U10A

#

$spe=$ARGV[3];
if($spe eq "fly")
{
	open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
	while(<IN>)
	{
	   if (/>(.+) type/)
	   {
	      $chr="chr$1";
	   }
	   else
	   {
	      chomp;
	      $genome{$chr}=$_;
	   }
	}
}
elsif($spe eq "bombyx")
{
	$fastafile="/home/wangw1/pipeline_bm/common/silkgenome.formatted.fa";
	open IN, $fastafile or die "Fail to open $fastafile: $!";
	while(<IN>)
	{
	   if (/>(.+)\s*\//) #this is specific for the case: >nscaf100 /length=4083 /lengthwogaps=4073
	   {
	      $chr="$1";
	      @c=split(/ /,$chr);
	   }
	   else
	   {
	      chomp;
	      $genome{$c[0]}=$_;
	   }
	}
}



my @matchedpairs=("AT","TA","GC","CG");
my @unmatchedpairs=("AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");
my @pairs=("AT","TA","GC","CG","AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");

my $OUTDIR=$ARGV[4];
my $indexFlag=$ARGV[5];

my %targetpf=();
my %guidepf=(); 

my %targetpfsplit=();
my %guidepfsplit=();

my %total=();

if($indexFlag)
{
	# make bowtie index
	for ($i=0; $i<$ARGV[2]; $i++)
	{
		my $file=fileparse($ARGV[$i]);
	    @namefield=split(/\./,$file);
	    if($spe eq "fly")
		{$name=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
	    push @argos, $name;
	    $file=$name;
	   	my $gz = gzopen($ARGV[$i], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
	   	while($gz->gzreadline($_) > 0)
	   	{
			chomp;
			split(/\t/);
			#($chr,$start,$end,$strand,$seq,$reands,$ntm)= split(/\t/,$_);
			
			next if (length($_[4])>29 || length($_[4])<23);
			next if (/data/);
			
			$total{$file}+=$_[5]/$_[6]; #no nta in norm.bed format
			
			#store the seq only; for faster extract the reads number later
			$guidepf{$file}{substr($_[4],0,16)}+=$_[5]/$_[6];

			#store chr, 5'end and strand information separately for each guide 16nt prefix	      
			if($_[3] eq '+')
			{ 
				$fiveend=$_[1]-1; #0-based
				$guidepfsplit{$file}{substr($_[4],0,16)}{$_[0],$fiveend,$_[3]}+=$_[5]/$_[6]; #become 0-based from norm.bed format
		  	}
		  	else
		  	{
		  		$fiveend=$_[2]-1; #0-based
		  		$guidepfsplit{$file}{substr($_[4],0,16)}{$_[0],$fiveend,$_[3]}+=$_[5]/$_[6]; #become 0-based from norm.bed format
		  	}
	      
	      	for ($n=1;$n<=20;$n++)
	      	{
	        	if ($_[3] eq '+')
	         	{
		            $start=$_[1]+$n-17; # $_[1] is 1 based but $start is 0 based
		            $str=substr($genome{$_[0]},$start,16); #substr function is 0 based
		            $str=&revfa($str);
		            $targetpf{$file}{$n}{$str}+=$_[5]/$_[6];
		            
		            #store chr, 5'end and strand information separately for each target 16nt prefix
		            $fiveend=$start+16; #0-based 
		            $targetpfsplit{$file}{$n}{$str}{$_[0],$fiveend,"-"}+=$_[5]/$_[6]; #store the strand information for target strand
	        	}
	         	else
	         	{
		            $start=$_[2]-$n;
		            $str=substr($genome{$_[0]},$start,16);
		            $targetpf{$file}{$n}{$str}+=$_[5]/$_[6];
		            
		            #store chr, 5'end and strand information separately for each guide 16nt prefix	
		            $fiveend=$start; #0-based
		            $targetpfsplit{$file}{$n}{$str}{$_[0],$fiveend,"+"}+=$_[5]/$_[6];
	        	}#ifelse
	      	}#for
		} #while
		$gz->gzclose(); 
		if ($total{$file}>10)
		{
			$seqFile="$OUTDIR/$file.seq";
			if( ! -s $seqFile )
			{
				open OUT, ">$seqFile";
				foreach (keys %{$guidepf{$file}})
				{
					print OUT "$_\t$guidepf{$file}{$_}\n" if (length($_)==16);
				}
			}
			for ($n=1;$n<=20;$n++)
			{
				$fa="$OUTDIR/$file.ref.$n.fa";
				$indexb="$OUTDIR/$file.$n";
				open OUT, ">$fa";
				foreach (keys %{$targetpf{$file}{$n}})
				{
					print OUT ">$_\t$targetpf{$file}{$n}{$_}\n$_\n" if (length($_)==16);
				}
				`bowtie-build $fa $indexb && rm $fa`;
	
			}#for loop of the n
		}#if the total
	} #for loop of the file
} #if
else
{
	for ($i=0; $i<$ARGV[2]; $i++)
	{
   		my $file=fileparse($ARGV[$i]);
    	@namefield=split(/\./,$file);
    	if($spe eq "fly")
		{$name=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
    	push @argos, $name;
    	$file=$name;
   		my $gz = gzopen($ARGV[$i], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
   		while($gz->gzreadline($_) > 0)
   		{
			chomp;
			split(/\t/);  
			next if (length($_[4])>29 || length($_[4])<23);
			next if (/data/);
			$total{$file}+=$_[5]/$_[6];
			$guidepf{$file}{substr($_[4],0,16)}+=$_[5]/$_[6];
	      
			#store chr, 5'end and strand information separately for each guide 16nt prefix
			if($_[3] eq '+')
			{ 
		  		$fiveend=$_[1]-1;#0-based
	      		$guidepfsplit{$file}{substr($_[4],0,16)}{$_[0],$fiveend,$_[3]}+=$_[5]/$_[6]; #become 0-based from norm.bed format
			}
			else
			{
		  		$fiveend=$_[2]-1;#0-based
		  		$guidepfsplit{$file}{substr($_[4],0,16)}{$_[0],$fiveend,$_[3]}+=$_[5]/$_[6]; #become 0-based from norm.bed format
			}
	      
			for ($n=1;$n<=20;$n++)
			{
				if ($_[3] eq '+')
				{
		            $start=$_[1]+$n-17;
		            $str=substr($genome{$_[0]},$start,16);
		            $str=&revfa($str);
		            $targetpf{$file}{$n}{$str}+=$_[5]/$_[6];
		            
		            #store chr, 5'end and strand information separately for each target 16nt prefix
		            $fiveend=$start+16; #0-based
		            $targetpfsplit{$file}{$n}{$str}{$_[0],$fiveend,"-"}+=$_[5]/$_[6];
	            }
	         	else
	         	{
		            $start=$_[2]-$n;
		            $str=substr($genome{$_[0]},$start,16);
		            $targetpf{$file}{$n}{$str}+=$_[5]/$_[6];
		            
		            #store chr, 5'end and strand information separately for each target 16nt prefix
		            $fiveend=$start;#0-based
		            $targetpfsplit{$file}{$n}{$str}{$_[0],$fiveend,"+"}+=$_[5]/$_[6];
	            
	         	}#ifelse
			}#for loop of the n
		}#while
		$gz->gzclose(); 
	}#for loop of the file
}#else indexFlag

open PPZ, ">$OUTDIR/zscore.toofewreads.out";
# bowtie mapping and score calculating
for ($i=0; $i<$ARGV[2]; $i++)
{
	for ($j=0; $j<=$i; $j++)
	{
      
		$file1=fileparse($ARGV[$i]); 
		@namefield=split(/\./,$file1);
   	    if($spe eq "fly")
		{$name1=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name1=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
		$file1=$name1;
		$file2=fileparse($ARGV[$j]);
		@namefield=split(/\./,$file2);
   	    if($spe eq "fly")
		{$name2=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name2=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
		$file2=$name2;
		#modify the order of filename on 02-10-2014 to clearly indicate guide target   


   
		if ($total{$file1}<10 || $total{$file2}<10) { print "$file2-$file1\t-10\n"; print PPZ "$file2-$file1\t-10\n";} #only when consider per cluster or per transposon family
		else
		{
			&PPprocessing($file2,$file1);
			if($file1 ne $file2 ) #added on 11/14/2013
			{   				
			   &PPprocessing($file1,$file2);    
			}#when file1 and file2 are the same the second iteration is skipped
		}#ifelse: total reads>10 

	}#j
}#i
close(PPZ);



sub PPprocessing
{

				my ($guideStrandFile,	$targetStrandFile)=@_;		

				my $X=0; 
				my $Z=0; 
				my %score=();
				my $count_N=0;
				my %X0=(); 
				my %Z0=(); 
				my %allPairReads=(); 
				my %suv=(); 
				my %count_N0=();
				my %species=(); 
				my %speciesn10=();
				my %cisPairReads=(); my %cisPairReads=(); my %cisPairReads=(); my %cisPairReads=();my %transPair10Species=();
				
				
				open ZSCORE, ">$OUTDIR/$guideStrandFile.$targetStrandFile.zscore.out";
				open ZSCOREUA, ">$OUTDIR/$guideStrandFile.UA_VA.zscore.out";
				 
				open PPSCORE, ">$OUTDIR/$guideStrandFile.$targetStrandFile.pp";
				open PPSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.VA.pp";
				
				$ppseq1="$OUTDIR/$guideStrandFile.$targetStrandFile.ppseq";
				open PPSEQ, ">$OUTDIR/$guideStrandFile.$targetStrandFile.ppseq";
				
				foreach ($n=1;$n<=20;$n++)
				{
					# file1 as ref
					$indexb="$OUTDIR/$targetStrandFile.$n";
					$seqFile="$OUTDIR/$guideStrandFile.seq";
					$bowtieOut="$OUTDIR/$guideStrandFile.$targetStrandFile.$n.bowtie.out";
				 	`bowtie $indexb -r -a -v 1 -p 8 $seqFile --suppress 1,4,6,7 | grep + > $bowtieOut`;
				   	%NTM=();
				   	open IN, "$OUTDIR/$guideStrandFile.$targetStrandFile.$n.bowtie.out";
				   	while($line=<IN>)
				   	{
					   	chomp $line;
					   	@l=split(/\t/,$line);
					   	#($strand, $bowtieindex,$seq,$mm)=split(/\t/,$line); 
					    
					   	if ($l[3]=~/(\d+):/)
					   	{ next if ($1!=0);}  # no seed mismatches 0 based
					   	$NTM{$l[2]}++;
					   	#$count_N++;
				   	}
				   	close(IN);
				   	open IN, "$OUTDIR/$guideStrandFile.$targetStrandFile.$n.bowtie.out";
				   	while($line=<IN>)
				   	{
				      	chomp $line;
				      	@l=split(/\t/,$line);
				      	if ($l[3] eq "")
				      	{
				      	
				      	   $g_0_nt=substr($l[2],0,1); $t_9_nt=&revfa($g_0_nt);  ##here are different from pp8_q2_ww1.pl
					       #targetpf index; guidepf seq
					       
					       $guidetotal=$guidepf{$guideStrandFile}{$l[2]};
					       $targettotal=$targetpf{$targetStrandFile}{$n}{$l[1]};
					       
					       $allPairReads{$g_0_nt.$t_9_nt}{$n}+=$targettotal*$guidetotal/$NTM{$l[2]};
					       
					       $speciesn10{$g_0_nt.$t_9_nt}{$l[2]}=1 if ($n==10); ###
					       $species{$g_0_nt.$t_9_nt}{$n}{$l[2]}=1 ; #this was wrong, has to add {$n}, otherwise accumulative
					       #the sum of $cisPairSpecies and $transPairSpecies irrespective of coordinates
					       
					       $score{$n}+=$targettotal*$guidetotal/$NTM{$l[2]};
				      		
					       #separate trans-targets from cis-targets for perfect g1t10 pair  	
					       
					       #$guidepfsplit{$file}{substr($_[4],0,16)}{$_[0],$fiveend,$_[3]}+=$_[5]/$_[6]; #become 0-based from 1-based norm.bed format
					       #$targetpfsplit{$file}{$n}{$str}{$_[0],$fiveend,"+"}+=$_[5]/$_[6];#$fiveend is 0-based
					       
					       #how to sort guidepfsplit and targetpfsplit
					       #no need to sort if use hash key
					       
					       
					       foreach my $record (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}} )
					       {
					       		my ($chr,$gfiveend,$gstrand)=split($record);
					       		
					       		
					       		#instead of scanning each key of target strand, get the 9 nt distant target directly
					       		#if($gstrand eq "+")
					       		#{
					       		#	my $tfiveend=$gfiveend+9; #10 bp overlap but 9 nt distance
					       		#	my $tstrand="-";
					       		#}
					       		#else
					       		#{
					       		#	my $tfiveend=$gfiveend-9; #10 bp overlap but 9 nt distance
					       		#	my $tstrand="+";					       			
					       		#}
					       		###query sequence 16 nt prefix population including both guides and targets
					       		###As I stored the potential target information above, if the potential target exist,
					       		###then it should find the excat match(chr,fiveend,strand) in the query seq population
					       		
					       		###note: how to define cistargets more accurately?
					       		###in addition to excat match, what if just several nucleotides away?
					       		
					       		my $tfiveend=$gfiveend;
					       		my $tstrand=$gstrand;
					       		
					       		
					       		
					       		if($targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$chr,$tfiveend,$tstrand})
					       		{
					       			$cisPairReads{$g_0_nt.$t_9_nt}{$n}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$chr,$tfiveend,$tstrand}/$NTM{$l[2]};
					       			$transPairReads{$g_0_nt.$t_9_nt}{$n}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*($targettotal - $targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$chr,$tfiveend,$tstrand})/$NTM{$l[2]};
					       			
					       			$cisPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}=1 ; #cis pair species must only have one by coordinate definition
					       			
					       			#with the same guide 16 nt prefix,there might be multiple trans-targets with 16nt complementarity, (originally it was viewed only one)
					       			$transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}})-1 ;
					       			#trans PingPong pair in species
					       			$transPair10Species{$g_0_nt.$t_9_nt}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}})-1 if($n==10) ;
					       			
					       		}
					       		else
					       		{
					       			$transPairReads{$g_0_nt.$t_9_nt}{$n}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targettotal/$NTM{$l[2]};
					       			$transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}});
					       			#trans PingPong pair in species
					       			$transPair10Species{$g_0_nt.$t_9_nt}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}}) if($n==10) ;
					       		}
		
					       }
					       
					       
					       #store target seq from query populations, as it's from bowtie output, by default, it has a partner
					       print PPSEQ "$l[2]\n" if ($n==10);
				       
				       
				       
				      	}#perfect pair

				      	elsif ($l[3]=~/(\d+):(\w)>(\w)/)
				      	{
					       next if ($1!=0);  # allow 1mm at the 10th position of target strand
					       $t_9_nt=&revfa($2);
					       $allPairReads{$3.$t_9_nt}{$n}+=$targetpf{$targetStrandFile}{$n}{$l[1]}*$guidepf{$guideStrandFile}{$l[2]}/$NTM{$l[2]};
					       $speciesn10{$3.$t_9_nt}{$l[2]}=1 if ($n==10); ###species of seq pairs, not count different coordinates
					       $species{$3.$t_9_nt}{$n}{$l[2]}=1 ;
					       $score{$n}+=$targetpf{$targetStrandFile}{$n}{$l[1]}*$guidepf{$guideStrandFile}{$l[2]}/$NTM{$l[2]}; ###species of pairs taking the varied coordinates
					       print PPSEQ "$l[2]\n" if ($n==10);
					       
					       
					       
					       $transPairReads{$3.$t_9_nt}{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targettotal/$NTM{$l[2]};
					       
					       $transPairSpecies{$3.$t_9_nt}{$n}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}});
					       
					       $transPair10Species{$3.$t_9_nt}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}}) if($n==10);
				      	}
				   } #while
				   close(IN);
				   
				   #total Ping-Pong
				   
				   	 $score{$n}=0 if (!exists $score{$n});
				     print PPSCORE "$n\t$score{$n}\n";
				     $count_N++ if ($score{$n}>0);
				   
				   #Ping-Pong score according to different G1T10 pairs
				   #for matched pairs, cis and trans  
				     foreach my $p (@matchedpairs)
				     {
					     $allPairReads{$p}{$n}=0 if (!exists $allPairReads{$p}{$n});
					     $n_of_species=scalar (keys %{$species{$p}{$n}});
					     
					     
					     my $n_of_cisPairReads=0;
					     map {$n_of_cisPairReads+=$_} values %{$cisPairReads{$p}{$n}} ;
					     my $n_of_cisPairSpecies=0;
					     $n_of_cisPairSpecies=scalar (keys %{$cisPairSpecies{$p}{$n}});					    
					     my $n_of_cisPairSpecies_cor=0;
					     map {$n_of_cisPairSpecies_cor+=$_} values %{$cisPairSpecies{$p}{$n}} ;
					     
					     
					     my $n_of_transPairReads=0;
					     map {$n_of_transPairReads+=$_} values %{$transPairReads{$p}{$n}} ;
					     my $n_of_transPairSpecies=0;					     					   
					     $n_of_transPairSpecies=scalar (keys %{$transPairSpecies{$p}{$n}});
					     my $n_of_transPairSpecies_cor=0;
					     map {$n_of_transPairSpecies_cor+=$_} values %{$transPairSpecies{$p}{$n}};
					     
					     print PPSCOREUA "$n\tcis\t$p\t$n_of_cisPairSpecies_cor\t$n_of_cisPairSpecies\t$n_of_cisPairReads\t$n_of_species\t$allPairReads{$p}{$n}\n";
					     print PPSCOREUA "$n\ttrans\t$p\t$n_of_transPairSpecies_cor\t$n_of_transPairSpecies\t$n_of_transPairReads\t$n_of_species\t$allPairReads{$p}{$n}\n";
					     
					     $count_N0{$p}++ if ($transPairSpecies{$p}{$n}>0);
				     }
				     #for unmatched pairs, only trans
				     foreach my $p (@unmatchedpairs)
				     {
					     $allPairReads{$p}{$n}=0 if (!exists $allPairReads{$p}{$n});
					     $n_of_species=scalar (keys %{$species{$p}{$n}});
					     
					     my $n_of_transPairReads=0;
					     map {$n_of_transPairReads+=$_} values %{$transPairReads{$p}{$n}} ;					     					   
					     $n_of_transPairSpecies=scalar (keys %{$transPairSpecies{$p}{$n}});
					     my $n_of_transPairSpecies_cor=0;
					     map {$n_of_transPairSpecies_cor+=$_} values %{$transPairSpecies{$p}{$n}};
					     
					     print PPSCOREUA "$n\ttrans\t$p\t$n_of_transPairSpecies_cor\t$n_of_transPairSpecies\t$n_of_transPairReads\t$n_of_species\t$allPairReads{$p}{$n}\n";
					     
					     $count_N0{$p}++ if ($transPairSpecies{$p}{$n}>0);
				     }
				     
				     
				     
				     
				}#n=1..20
				   $X=$score{10}; delete $score{10};
				   $std=&std(values %score); 
				   if ($std>0 && $count_N>=5) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
				   print ZSCORE "$guideStrandFile\-$targetStrandFile\t$Z\t";
				   print "$guideStrandFile\-$targetStrandFile\t$Z\t";
				   
				   #Z-score for all transpairs
				   foreach my $p (@pairs)
				   {
				    #$std=0;
				    $X0{$p}=$transPairSpecies{$p}{10}; delete $transPairSpecies{$p}{10};
				    $std0=&std(values %{$transPairSpecies{$p}});
				    if ($std0>0 && $count_N0{$p}>=5) { $Z0{$p}=($X0{$p}-&mean(values %{$transPairSpecies{$p}}))/$std0;} else {$Z0{$p}=-10;}
				    $n_of_species=scalar (keys %{$transPair10Species{$p}}); #total number of species irrespective of coordinates
				    my $n_of_species_cor=0;
				    map {$n_of_species_cor+=$_} values %{$transPair10Species{$p}};#total number of species respective of coordinates
				    #how to normalize $X0{$p}?
				    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\t$p\t$Z0{$p}\t$X0{$p}\t$n_of_species\t$n_of_species_cor\n"; ##file2 is the guide and file1 is the target
				   }
				   
				   $ppseq="$OUTDIR/$guideStrandFile.$targetStrandFile.ppseq";
				   $seqFile="$OUTDIR/$guideStrandFile.seq";
				   $NofPPreads=`match.pl $ppseq $seqFile | sumcol+ 2`; chomp($NofPPreads);
				   
				   if ($Z!=-10) 
				   {
				   	print ZSCORE "$X\t$NofPPreads\t$total{$guideStrandFile}",$NofPPreads/$total{$guideStrandFile},"\t",$X*1000000000000/$total{$targetStrandFile}/$total{$guideStrandFile};
				   }
				   else
				   {
				   	print ZSCORE "NA\tNA\tNA\tNA\tNA\n";
				   }
				   
				close(PPSEQ);
				close(PPSCOREUA);
				close(PPSCORE);
				close(ZSCOREUA);
				close(ZSCORE);  	
}

#remove intermediate files
#`rm *.ebwt`;
#`rm *.bowtie.out`;
