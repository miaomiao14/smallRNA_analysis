#!/usr/bin/perl

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "restrict_digts.pm";
require "Jia.pm";
use File::Basename;
use Compress::Zlib;
use List::Util qw(sum);
# rule is p1 and p17-21 doesn't need to pair but p2-16 need
# simplifized version with prefix 16nt
# input as norm.bed ( no header)

#this version is staring form the target, and look for guides (converting the target information to guides), make index of the potential guides

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

#02/18/2014
#add bed format

#02/23/2014
#normalize species, different from 02102014 version

#02/25/2014
#fix a bug for base fraction; target strand


#03/02/2014
#add winsize, prefix as parameters

#03/04
#add bowtie index dir

if(scalar(@ARGV)<6)
{
        usage();
}

my $parameters={};#initialize a reference of a hash

&parse_command_line($parameters, @ARGV);

my $inputfile1=fileparse($parameters->{in_file1});
my @inputfiles=();
push @inputfiles, $parameters->{in_file1};
my $inputfile2="";
if($parameters->{in_file2})
{
	$inputfile2=fileparse($parameters->{in_file2});
	push @inputfiles, $parameters->{in_file2};	
}
my $numOfInput=$parameters->{num};
my $spe=$parameters->{species};
my $OUTDIR=$parameters->{outdir};
my $BOUTDIR=$parameters->{indexoutdir};
my $MOUTDIR=$parameters->{mappingoutdir};
my $QOUTDIR=$parameters->{queryseqoutdir};
my $indexFlag=$parameters->{indexflag};
my $fileFormat=$parameters->{format};
my $wsize=$parameters->{winsize};
my $basep=$parameters->{complementarity};
my $fastafile=$parameters->{fa};

if($spe eq "fly")
{
	open IN, $fastafile;
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
elsif($spe eq "mouse" or "gfp")
{
	open IN, $fastafile or die "Fail to open $fastafile: $!";
	while(my $line=<IN>)
	{
	   if ($line=~/>(.+)/) #>chr1
	   {
	      $chr="$1";
	   }
	   else
	   {
	      chomp $line;
	      $genome{$chr}=$line;
	   }
	}
	close(IN);
}


my @matchedpairs=("AT","TA","GC","CG");
my @unmatchedpairs=("AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");
my @pairs=("AT","TA","GC","CG","AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");



my %targetpf=();
my %guidepf=(); 

my %targetpfsplit=();
my %guidepfsplit=();

my %total=();

my %totalFirstBase=();
my %totalTenthBase=();
#my %totalExpectedTenthBase=();

my %totalTenthBaseTrial=();
my %totalFirstBaseSpecies=();
my %totalTenthBaseSpecies=();

#main, preprocessing
if($indexFlag)
{
	# make bowtie index
	for ($i=0; $i<$numOfInput; $i++)
	{
		my $file=fileparse($inputfiles[$i]);
	    @namefield=split(/\./,$file);
	    if($spe eq "fly" or "gfp")
		{$name=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
	    if($spe eq "mouse")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13]."_".$namefield[6];}
	    push @argos, $name;
	    $file=$name;
	    
	    &InputFileProcessing($inputfiles[$i],$file);
	     
		if ($total{$file}>10)
		{
			$seqFile="$QOUTDIR/$file.seq";
			if( ! -s $seqFile )#test the existence of file
			{
				open OUT, ">$seqFile";
				foreach my $prefix (keys %{$guidepf{$file}})
				{
					print OUT "$prefix\t$guidepf{$file}{$prefix}\n" if (length($prefix)==$basep);
				}
			}
			for ($n=0;$n<$wsize;$n++)
			{
				$fa="$OUTDIR/$file.ref.$n.fa";
				$indexb="$BOUTDIR/$file.$basep.$n";
				if(! -s $indexb) #test the existence of file
				{
					open OUT, ">$fa";
					foreach my $prefix (keys %{$targetpf{$file}{$n}})
					{
						print OUT ">$prefix\t$targetpf{$file}{$n}{$prefix}\n$prefix\n" if (length($prefix)==$basep);
					}
					`bowtie-build $fa $indexb && rm $fa`;
				}
	
			}#for loop of the n
		}#if the total
	} #for loop of the file
} #if
else #if indexFlag
{
	for ($i=0; $i<$numOfInput; $i++)
	{
   		my $file=fileparse($inputfiles[$i]);
    	@namefield=split(/\./,$file);
    	if($spe eq "fly" or "gfp")
		{$name=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
	    if($spe eq "mouse")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13]."_".$namefield[6];}
    	push @argos, $name;
    	$file=$name;
   		&InputFileProcessing($inputfiles[$i],$file);
	}#for loop of the file
}#else indexFlag

#main
open PPZ, ">$OUTDIR/zscore.toofewreads.out";
# bowtie mapping and score calculating
for ($i=0; $i<$numOfInput; $i++)
{
	for ($j=0; $j<=$i; $j++)
	{
      
		      
		$file1=fileparse($inputfiles[$i]); 
		@namefield=split(/\./,$file1);
   	    if($spe eq "fly" or "gfp")
		{$name1=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name1=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
	    if($spe eq "mouse")
	    {$name1=$namefield[2]."_".$namefield[12]."_".$namefield[13]."_".$namefield[6];}
		$file1=$name1;
		$file2=fileparse($inputfiles[$j]);
		@namefield=split(/\./,$file2);
   	    if($spe eq "fly" or "gfp")
		{$name2=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name2=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
	    if($spe eq "mouse")
	    {$name2=$namefield[2]."_".$namefield[12]."_".$namefield[13]."_".$namefield[6];}
		$file2=$name2;
		#modify the order of filename on 02-10-2014 to clearly indicate guide target   


   
		if ($total{$file1}<10 || $total{$file2}<10) {print PPZ "$file2-$file1\t-10\n";} #only when consider per cluster or per transposon family
		else
		{
			&PingPongProcessing($file2,$file1);
			if($file1 ne $file2 ) #added on 11/14/2013
			{   				
			   &PingPongProcessing($file1,$file2);    
			}#when file1 and file2 are the same the second iteration is skipped
		}#ifelse: total reads>10 

	}#j
}#i
close(PPZ);

sub InputFileProcessing
{
	my ($inputfile,$file)= @_;
	my $gz = gzopen($inputfile, "rb") or die "Cannot open $inputfile: $gzerrno\n" ;
   	while($gz->gzreadline($line) > 0)
   	{
		chomp $line;
		#split(/\t/);
		my $chr;
		my $bedstart;
		my $bedend;
		my $strand;
		my $seq;
		my $reads;
		my $ntm;
		if($fileFormat eq "normbed")
		{ 
			($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$line);
		}
		if($fileFormat eq "bed")
		{
			($chr,$bedstart,$bedend,$seq,$ntmreads,$strand)= split(/\t/,$line);
			$reads=$ntmreads;
			$ntm=1;
		}
		
		if($inputfile=~/SRA/){next if (length($seq)>29 || length($seq)<23);}
		next if (/data/);
		
		$total{$file}+=$reads/$ntm; #no nta in norm.bed format
		#store chr, 5'end and strand information separately for each query 20nt prefix, the populations to find the guide
		my $fiveend=0;
		my $dnaseq="";	      
		if($strand eq '+')
		{	
			
			if($fileFormat eq "bed")
			{ 
				$fiveend=$bedstart; #0-based,closed
			}
			if($fileFormat eq "normbed")
			{
				$fiveend=$bedstart-1;#convert to 0-based,closed

			}
			$dnaseq=substr($genome{$chr},$fiveend,$basep);
			$totalFirstBase{$file}{substr($dnaseq,0,1)}{$dnaseq}+=$reads/$ntm;
		
			#store the seq of guide 20nt prefix only; for faster extract the reads number later
			$guidepf{$file}{$dnaseq}+=$reads/$ntm;
			$guidepfsplit{$file}{substr($dnaseq,0,$basep)}{"$chr,$fiveend,$strand"}+=$reads/$ntm; #become 0-based from norm.bed format
			

	  	}
	  	else
	  	{
	  		if($fileFormat eq "bed")
			{
	  			$fiveend=$bedend; #open
			}
			if($fileFormat eq "normbed")
			{
				$fiveend=$bedend;#convert to bedformat,open
				
				
			}
			my $seqstart=$fiveend-$basep;
     		my $seqtemp=substr($genome{$chr},$seqstart,$basep);
     		$dnaseq=&revfa($seqtemp);
				
			$totalFirstBase{$file}{substr($dnaseq,0,1)}{$dnaseq}+=$reads/$ntm;
		
				#store the seq of guide 20nt prefix only; for faster extract the reads number later
			$guidepf{$file}{$dnaseq}+=$reads/$ntm;
			$guidepfsplit{$file}{substr($dnaseq,0,$basep)}{"$chr,$fiveend,$strand"}+=$reads/$ntm; #become 0-based from norm.bed format

	  	}
      	for (my $n=0;$n<$wsize;$n++)
      	{
      		my $start=0;
      		my $fiveend=0;
        	if ($strand eq '+') #target strand information
         	{
         		
         		if($fileFormat eq "bed")
         		{
	            	$start=$bedstart+$n+1-$basep; # $bedstart is 0 based and $start is 0 based, the intermediate end $bedstart+$n+1 is open
	            	$fiveend=$start+$basep; #open
         		}
         		if($fileFormat eq "normbed")
         		{
         			$start=$bedstart+$n-$basep;
         			$fiveend=$start+$basep;#convert to bed, open
         		}
	            my $str=substr($genome{$chr},$start,$basep); #substr function is 0 based
	            $str=&revfa($str);
	            $targetpf{$file}{$n}{$str}+=$reads/$ntm;
	            
	            #store chr, 5'end and strand information separately for each guide 16nt prefix
	            my $tstrand="-";
	            $targetpfsplit{$file}{$n}{$str}{"$chr,$fiveend,$tstrand"}+=$reads/$ntm; #store the strand information for guide strand
	            #my $indexStart=$start;
	            
	            $totalTenthBaseTrial{$file}{$n}{substr($dnaseq,$n,1)}{$str}{substr($dnaseq,0,$basep)}+=$reads/$ntm;
	            
        	}
         	else
         	{
         		if($fileFormat eq "bed")
         		{
		            $start=$bedend-$n-1; #closed
		            $fiveend=$start; #0-based
         		}
         		if($fileFormat eq "normbed")
         		{
         			$start=$bedend-$n-1; #closed
         			$fiveend=$start;	
         		}
	            my $str=substr($genome{$chr},$start,$basep);
	            $targetpf{$file}{$n}{$str}+=$reads/$ntm;
	            my $tstrand="+";
	            #store chr, 5'end and strand information separately for each guide 16nt prefix	
	            $targetpfsplit{$file}{$n}{$str}{"$chr,$fiveend,$tstrand"}+=$reads/$ntm;
	            
	            $totalTenthBaseTrial{$file}{$n}{substr($dnaseq,$n,1)}{$str}{substr($dnaseq,0,$basep)}+=$reads/$ntm;
	            
        	}#ifelse
      	}#for
	} #while
	$gz->gzclose();
}

sub PingPongProcessing
{

	my ($guideStrandFile,	$targetStrandFile)=@_;		

	my $X=0; 
	my $Z=0; 
	my %score=();
	my $count_N=0;
	my %allPairReads=(); 
	my %count_N0=();
	my %species=(); 
	my %speciesn10=();
	
	my %cisPairSpecies=();
	my %cisPair10Species=();
	my %cisPairReads=();
	
	my %transPairSpecies=();
	my %transPair10Species=();
	my %transPairReads=();
	
	
	my %pp6cisPairSpecies=();
	my %pp6cisPair10Species=();
	my %pp6cisPairReads=();
	
	my %pp8allPairSpecies=();
	my %pp8allPair10Species=();
	my %pp8allPairReads=();
	
	my %transallPairSpecies=();
	my %transallPair10Species=();
	my %transallPairReads=();
	
	

	
	open ZSCORE, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.zscore.out";
	open ZSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.zscore.out";
	 
	open PPSCORE, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.pp";
	open PPSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.pp";
	
	open PPUAFRACTION, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.base.fraction.txt";
	
	open PPGSEQ, ">$OUTDIR/$guideStrandFile.$basep.prefix.G.ppseq";
	open PPTSEQ, ">$OUTDIR/$targetStrandFile.$basep.prefix.T.ppseq";
	
	foreach ($n=0;$n<$wsize;$n++)
	{
		my %pairedFirstBase=();
		my %pairedTenthBase=();
		# file1 as ref
		$indexb="$BOUTDIR/$targetStrandFile.$basep.$n";
		$seqFile="$QOUTDIR/$guideStrandFile.seq";
		$bowtieOut="$MOUTDIR/$guideStrandFile.$targetStrandFile.$basep.$n.bowtie.out";
	 	`[ ! -f $bowtieOut ] && bowtie $indexb -r -a -v 1 -p 8 $seqFile --suppress 1,4,6,7 | grep + > $bowtieOut`;
	   	my %NTM=();
	   	open IN, "$MOUTDIR/$guideStrandFile.$targetStrandFile.$basep.$n.bowtie.out";
	   	while(my $line=<IN>)
	   	{
		   	chomp $line;
		   	@l=split(/\t/,$line);
		   	#($strand, $bowtieindex,$seq,$mm)=split(/\t/,$line); 
		    
		   	if ($l[3]=~/(\d+):/)
		   	{ next if ($1!=0);}  # no seed mismatches 0 based
		   	$NTM{$l[2]}++;
	   	}
	   	close(IN);
	   	open IN, "$MOUTDIR/$guideStrandFile.$targetStrandFile.$basep.$n.bowtie.out";
	   	while(my $line=<IN>)
	   	{
	      	chomp $line;
	      	@l=split(/\t/,$line);
	      	
	      	my $guidetotal=$guidepf{$guideStrandFile}{$l[2]};
		    my $targettotal=$targetpf{$targetStrandFile}{$n}{$l[1]};
	      	
	      	my $gttotal=$guidetotal*$targettotal;
	      	
	      	my $nGcor=scalar (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}});
		    my $nTcor=scalar (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}});
		    
		    my $nnGcorTcor=$nGcor*$nTcor;
	      	
	      	if ($l[3] eq "")
	      	{
	      	
	      	   $g_0_nt=substr($l[2],0,1); $t_9_nt=&revfa($g_0_nt);  ##here are different from pp8_q2_ww1.pl
		       #targetpf index; guidepf seq
			   $score{$n}+=$gttotal/$NTM{$l[2]}; #total pp8 ppscore
			   
			   #how many of species start with U?
  
		       	$pairedFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]}=$totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]} if($totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]});
			   	# expect to give the same results as %pairedTenthBase
				
				#$pairedExpectedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$l[1]}=$totalExpectedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$l[1]} if($totalExpectedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$l[1]});#should not be accumulative			 
				#compare to the results by ppseq files, the key(l[1],expected guides) reduce the number of targets from different locations (not every targets have the same expected guides) 
				#however, compare to the %totalTenthBaseTrial, it has more reads, it has all targets link to the expected guides	
				if($totalTenthBaseTrial{$targetStrandFile}{$n}{$t_9_nt}{$l[1]})
				{
					foreach my $seq (keys %{$totalTenthBaseTrial{$targetStrandFile}{$n}{$t_9_nt}{$l[1]}})
					{
						$pairedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$seq}+=$totalTenthBaseTrial{$targetStrandFile}{$n}{$t_9_nt}{$l[1]}{$seq};
						print PPTSEQ "$seq\n" if ($n==9);
					}
				}
				print PPGSEQ "$l[2]\n" if ($n==9);
		        
			   $species{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]} ; #this was wrong, has to add {$n}, otherwise accumulative
		       #the sum of $cisPairSpecies and $transPairSpecies taking coordinates
		       $speciesn10{$g_0_nt.$t_9_nt}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]} if ($n==9); ###
		       
		       $allPairReads{$g_0_nt.$t_9_nt}{$n}+=$gttotal/$NTM{$l[2]};
		       
		       $pp8allPairReads{$n}{$l[2]}+=$gttotal/$NTM{$l[2]}; 
		       $pp8allPairSpecies{$n}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]};#no g1t10 as key compare to the species hash, that's why fewer species from this pp8allPairSpecies
		       $pp8allPair10Species{$l[2]}+=$nnGcorTcor/$NTM{$l[2]} if($n==9) ;
		       		     		       
		       foreach my $record (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}} )
		       {
		       		my ($chr,$gfiveend,$gstrand)=split(/,/,$record);

		       		###note: how to define cistargets more accurately?
		       		###in addition to excat match, what if just several nucleotides away?
		       		
		       		my $tfiveend=$gfiveend;
		       		my $tstrand=$gstrand;

		       		if($targetpfsplit{$targetStrandFile}{$n}{$l[1]}{"$chr,$tfiveend,$tstrand"})
		       		{
		       			$cisPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{"$chr,$tfiveend,$tstrand"}/$NTM{$l[2]};
		       			
		       			$cisPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=1/$NTM{$l[2]} ; #cis pair species must only have one by coordinate definition; but for different record, it has different cis pair
		       			$cisPair10Species{$g_0_nt.$t_9_nt}{$l[2]}+=1/$NTM{$l[2]} if($n==9) ;
		       			
		       			$pp6cisPairReads{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{"$chr,$tfiveend,$tstrand"};#for pp6 does not care mismatched pairs, no need /$NTM{$l[2]}
		       			$pp6cisPairSpecies{$n}{$l[2]}+=1/$NTM{$l[2]};
		       			$pp6cisPair10Species{$l[2]}+=1/$NTM{$l[2]} if($n==9) ;
		       			
		       			#with the same guide 16 nt prefix,there might be multiple trans-targets with 16nt complementarity, (originally it was viewed only one)
		       			#trans PingPong pair in species
		       			$transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=($nTcor-1)/$NTM{$l[2]} ;
		       			$transPair10Species{$g_0_nt.$t_9_nt}{$l[2]}+=($nTcor-1)/$NTM{$l[2]} if($n==9) ;
		       			#trans PingPong pair in reads
		       			$transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*($targettotal - $targetpfsplit{$targetStrandFile}{$n}{$l[1]}{"$chr,$tfiveend,$tstrand"})/$NTM{$l[2]};
		       			
		       			$transallPairSpecies{$n}{$l[2]}+=($nTcor-1)/$NTM{$l[2]} ;
						$transallPair10Species{$n}{$l[2]}+=($nTcor-1)/$NTM{$l[2]} if($n==9);
						$transallPairReads{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*($targettotal - $targetpfsplit{$targetStrandFile}{$n}{$l[1]}{"$chr,$tfiveend,$tstrand"})/$NTM{$l[2]};
		       			
		       			
		       		}
		       		else
		       		{
		       			#trans PingPong pair in species
		       			$transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$nTcor/$NTM{$l[2]};
		       			#trans PingPong pair in reads
		       			$transPair10Species{$g_0_nt.$t_9_nt}{$l[2]}+=$nTcor/$NTM{$l[2]} if($n==9) ;
		       			$transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targettotal/$NTM{$l[2]};
		       			
		       			$transallPairSpecies{$n}{$l[2]}+=$nTcor/$NTM{$l[2]};
						$transallPair10Species{$n}{$l[2]}+=$nTcor/$NTM{$l[2]} if($n==9);
						$transallPairReads{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targettotal/$NTM{$l[2]};
		       		}

		       }
		       					       
		       #store target seq from query populations, as it's from bowtie output, by default, it has a partner

	      	}#perfect pair

	      	elsif ($l[3]=~/(\d+):(\w)>(\w)/)
	      	{

		       next if ($1!=0);  # allow 1mm at the 10th position of target strand
		       
		       $g_0_nt=$3;
		       $t_9_nt=&revfa($2);
		       
		       $score{$n}+=$gttotal/$NTM{$l[2]}; #total pp8 ppscore
		       
		         
		       	$pairedFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]}=$totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]} if($totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]});
			   	# expect to give the same results as %pairedTenthBase
				
				#$pairedExpectedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$l[1]}=$totalExpectedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$l[1]} if($totalExpectedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$l[1]});#should not be accumulative			 
				#compare to the results by ppseq files, the key(l[1],expected guides) reduce the number of targets from different locations (not every targets have the same expected guides) 
				#however, compare to the %totalTenthBaseTrial, it has more reads, it has all targets link to the expected guides	
				if($totalTenthBaseTrial{$targetStrandFile}{$n}{$t_9_nt}{$l[1]})
				{
					foreach my $seq (keys %{$totalTenthBaseTrial{$targetStrandFile}{$n}{$t_9_nt}{$l[1]}})
					{
						$pairedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$seq}+=$totalTenthBaseTrial{$targetStrandFile}{$n}{$t_9_nt}{$l[1]}{$seq};
						print PPTSEQ "$seq\n" if ($n==9);
					}
				}
				print PPGSEQ "$l[2]\n" if ($n==9);

		        		       
		       $species{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]}  ;###species of seq pairs, not count different coordinates
		       $speciesn10{$g_0_nt.$t_9_nt}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]}  if ($n==9); ###species of seq pairs, not count different coordinates
		       
		       $allPairReads{$g_0_nt.$t_9_nt}{$n}+=$gttotal/$NTM{$l[2]};###reads of seq pairs, not count different coordinates
				
			   
		       	
		       #trans PingPong pair in species
		       $transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]};
		       $transPair10Species{$g_0_nt.$t_9_nt}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]} if($n==9);
		       #trans PingPong pair in reads
		       $transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$gttotal/$NTM{$l[2]};
		       
		       $pp8allPairSpecies{$n}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]};
		       $pp8allPair10Species{$l[2]}+=$nnGcorTcor/$NTM{$l[2]} if($n==9) ;
		       $pp8allPairReads{$n}{$l[2]}+=$gttotal/$NTM{$l[2]};
		       
		       $transallPairSpecies{$n}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]};
			   $transallPair10Species{$n}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]} if($n==9);
			   $transallPairReads{$n}{$l[2]}+=$gttotal/$NTM{$l[2]};

			}
		} #while
		close(IN);
		
		 #$firstBaseFraction

	   my %pairedFirstBaseReads=();
	   my $pairedFirstBaseReadsTotal=0;
	   my %pairedFirstBaseReadsF=();
	   my %pairedFirstBaseSpecies=();
	   my $pairedFirstBaseSpeciesTotal=0;
	   my %pairedFirstBaseSpeciesF=();
	   
	   my %totalFirstBaseReads=();
	   my $totalFirstBaseReadsTotal=0;
	   my %totalFirstBaseReadsF=();
	   my %totalFirstBaseSpecies=();
	   my $totalFirstBaseSpeciesTotal=0;
	   my %totalFirstBaseSpeciesF=();
	   
	   foreach my $b (keys %{$pairedFirstBase{$guideStrandFile}})
	   {
	   		map {$pairedFirstBaseReads{$b}+=$_} values  %{$pairedFirstBase{$guideStrandFile}{$b}};
	   		$pairedFirstBaseReadsTotal+=$pairedFirstBaseReads{$b};
	   		$pairedFirstBaseSpecies{$b}=scalar (keys  %{$pairedFirstBase{$guideStrandFile}{$b}});
	   		$pairedFirstBaseSpeciesTotal+=$pairedFirstBaseSpecies{$b};
	   		
	   		#$totalFirstBase{$file}{substr($seq,0,1)}{substr($seq,0,20)}+=$reads/$ntm;
	   		#total
	   		map {$totalFirstBaseReads{$b}+=$_} values %{$totalFirstBase{$guideStrandFile}{$b}};
	   		$totalFirstBaseReadsTotal+= $totalFirstBaseReads{$b};
	   		$totalFirstBaseSpecies{$b}=scalar (keys %{$totalFirstBase{$guideStrandFile}{$b}});
	   		$totalFirstBaseSpeciesTotal+= $totalFirstBaseSpecies{$b};
	   }
	   foreach my $b (keys  %pairedFirstBaseReads)
	   {			
	   		$pairedFirstBaseReadsF{$b}=$pairedFirstBaseReads{$b}/$pairedFirstBaseReadsTotal;
	   		$pairedFirstBaseReadsF{$b}=&restrict_num_decimal_digits($pairedFirstBaseReadsF{$b},4);
	   		$pairedFirstBaseSpeciesF{$b}=$pairedFirstBaseSpecies{$b}/$pairedFirstBaseSpeciesTotal;
	   		$pairedFirstBaseSpeciesF{$b}=&restrict_num_decimal_digits($pairedFirstBaseSpeciesF{$b},4);
	   		
	   		$totalFirstBaseReadsF{$b}=$totalFirstBaseReads{$b}/$totalFirstBaseReadsTotal;
	   		$totalFirstBaseReadsF{$b}=&restrict_num_decimal_digits($totalFirstBaseReadsF{$b},4);
	   		$totalFirstBaseSpeciesF{$b}=$totalFirstBaseSpecies{$b}/$totalFirstBaseSpeciesTotal;
	   		$totalFirstBaseSpeciesF{$b}=&restrict_num_decimal_digits($totalFirstBaseSpeciesF{$b},4);
	   }
	   
	   	$pairedFirstBaseSpeciesTotal=&restrict_num_decimal_digits($pairedFirstBaseSpeciesTotal,4);
		$pairedFirstBaseReadsTotal=&restrict_num_decimal_digits($pairedFirstBaseReadsTotal,4);
	   	$totalFirstBaseSpeciesTotal=&restrict_num_decimal_digits($totalFirstBaseSpeciesTotal,4);
		$totalFirstBaseReadsTotal=&restrict_num_decimal_digits($totalFirstBaseReadsTotal,4);
	   #tenthBaseFraction
	   my %pairedTenthBaseReads=();
	   my $pairedTenthBaseReadsTotal=0;
	   my %pairedTenthBaseReadsF=();
	   my %pairedTenthBaseSpecies=();
	   my $pairedTenthBaseSpeciesTotal=0;
	   my %pairedTenthBaseSpeciesF=();
	   
	   
	   my %totalTenthBaseReads=();
	   my $totalTenthBaseReadsTotal=0;
	   my %totalTenthBaseReadsF=();
	   my %totalTenthBaseSpecies=();
	   my $totalTenthBaseSpeciesTotal=0;
	   my %totalTenthBaseSpeciesF=();
	   
	   foreach my $b (keys %{$pairedTenthBase{$targetStrandFile}{$n}})
	   {
	   		map {$pairedTenthBaseReads{$b}+=$_} values  %{$pairedTenthBase{$targetStrandFile}{$n}{$b}};
	   		$pairedTenthBaseReadsTotal+=$pairedTenthBaseReads{$b};
	   		#species
	   		$pairedTenthBaseSpecies{$b}=scalar (keys  %{$pairedTenthBase{$targetStrandFile}{$n}{$b}});
	   		$pairedTenthBaseSpeciesTotal+=$pairedTenthBaseSpecies{$b};
	   		
	   		
	   }

	  	   
	   my %totalTenthBaseTemp=();
	   foreach my $b (keys %{$totalTenthBaseTrial{$targetStrandFile}{$n}})
	   {
	   		
	   		foreach my $expTargets (keys %{$totalTenthBaseTrial{$targetStrandFile}{$n}{$b}})
	   		{
		   			#$totalTenthBase{$file}{substr($seq,0,1)}{substr($seq,0,20)}+=$reads/$ntm;
		   		foreach my $realTargets (keys %{$totalTenthBaseTrial{$targetStrandFile}{$n}{$b}{$expTargets}})
		   		{
		   		$totalTenthBaseTemp{$targetStrandFile}{$n}{$b}{$realTargets}+=$totalTenthBaseTrial{$targetStrandFile}{$n}{$b}{$expTargets}{$realTargets};
		   		}

		   		
	   		}
	   		#reads
		   	map {$totalTenthBaseReads{$b}+=$_} values %{$totalTenthBaseTemp{$targetStrandFile}{$n}{$b}};
		   	$totalTenthBaseReadsTotal+= $totalTenthBaseReads{$b};
		   	#species
		   	$totalTenthBaseSpecies{$b}+=scalar (keys %{$totalTenthBaseTemp{$targetStrandFile}{$n}{$b}});
	   		
	   		$totalTenthBaseSpeciesTotal+= $totalTenthBaseSpecies{$b};
	   }
	   

	   
	   #$pairedExpectedTenthBase
	   
	   foreach my $b (keys  %pairedTenthBaseReads)
	   {			
	   		$pairedTenthBaseReadsF{$b}=$pairedTenthBaseReads{$b}/$pairedTenthBaseReadsTotal;
	   		$pairedTenthBaseReadsF{$b}=&restrict_num_decimal_digits($pairedTenthBaseReadsF{$b},4);
	   		$pairedTenthBaseSpeciesF{$b}=$pairedTenthBaseSpecies{$b}/$pairedTenthBaseSpeciesTotal;
	   		$pairedTenthBaseSpeciesF{$b}=&restrict_num_decimal_digits($pairedTenthBaseSpeciesF{$b},4);
	   		
	   		$totalTenthBaseReadsF{$b}=$totalTenthBaseReads{$b}/$totalTenthBaseReadsTotal;
	   		$totalTenthBaseReadsF{$b}=&restrict_num_decimal_digits($totalTenthBaseReadsF{$b},4);
	   		$totalTenthBaseSpeciesF{$b}=$totalTenthBaseSpecies{$b}/$totalTenthBaseSpeciesTotal;
	   		$totalTenthBaseSpeciesF{$b}=&restrict_num_decimal_digits($totalTenthBaseSpeciesF{$b},4);
	   }

		$pairedTenthBaseSpeciesTotal=&restrict_num_decimal_digits($pairedTenthBaseSpeciesTotal,4);
		$pairedTenthBaseReadsTotal=&restrict_num_decimal_digits($pairedTenthBaseReadsTotal,4);
		

		$totalTenthBaseSpeciesTotal=&restrict_num_decimal_digits($totalTenthBaseSpeciesTotal,4);
		$totalTenthBaseReadsTotal=&restrict_num_decimal_digits($totalTenthBaseReadsTotal,4);
		
			  
		
		
		
		
		
		$m=$n+1;	
		my @bases=("A","C","G","T");
		print PPUAFRACTION "$m\ttotal\tg1\t$pairedFirstBaseSpeciesTotal\t$totalFirstBaseSpeciesTotal\t$pairedFirstBaseReadsTotal\t$totalFirstBaseReadsTotal\n";
		
		foreach my $b(@bases)
		{
			print PPUAFRACTION "$m\t$b\tg1\t$pairedFirstBaseSpeciesF{$b}\t$totalFirstBaseSpeciesF{$b}\t$pairedFirstBaseReadsF{$b}\t$totalFirstBaseReadsF{$b}\n";
			
		}
		print PPUAFRACTION "$m\ttotal\tt10\t$pairedTenthBaseSpeciesTotal\t$totalTenthBaseSpeciesTotal\t$pairedTenthBaseReadsTotal\t$totalTenthBaseReadsTotal\n";
		foreach my $b(@bases)
		{
			print PPUAFRACTION "$m\t$b\tt10\t$pairedTenthBaseSpeciesF{$b}\t$totalTenthBaseSpeciesF{$b}\t$pairedTenthBaseReadsF{$b}\t$totalTenthBaseReadsF{$b}\n";
		}

		if($m==10)
		{
			$ppgseq="$OUTDIR/$guideStrandFile.$basep.prefix.G.ppseq";#$basep.prefix.
			$seqFile="$OUTDIR/$guideStrandFile.seq";
			$ppgseqm="$OUTDIR/$guideStrandFile.$basep.prefix.G.ppseq.reads";
			`match.pl $ppgseq $seqFile >$ppgseqm`;
		
			my @bases=("A","C","G","T");
			#total
			my %totalG1guideStat=();
			my $totalG1guideStatRef=&g1Frac($seqFile);
			%totalG1guideStat=%{$totalG1guideStatRef};
			my $totalG1Species=0;
			my $totalG1Reads=0;
			my %totalG1SpeciesC=();
			my %totalG1ReadsC=();
			my %totalG1SpeciesF=();
			my %totalG1ReadsF=();
			foreach my $b (keys %totalG1guideStat)
			{
				$totalG1SpeciesC{$b}=scalar (keys %{$totalG1guideStat{$b}});
				$totalG1Species+=$totalG1SpeciesC{$b};
				map {$totalG1ReadsC{$b}+=$_} values %{$totalG1guideStat{$b}};
				$totalG1Reads+=$totalG1ReadsC{$b};
			}
		
			foreach my $b(keys %totalG1SpeciesC )
			{
				$totalG1SpeciesF{$b}=$totalG1SpeciesC{$b}/$totalG1Species;
				$totalG1SpeciesF{$b}=&restrict_num_decimal_digits($totalG1SpeciesF{$b},4);
				$totalG1ReadsF{$b}=$totalG1ReadsC{$b}/$totalG1Reads;
				$totalG1ReadsF{$b}=&restrict_num_decimal_digits($totalG1ReadsF{$b},4);
			}
			$totalG1Species=&restrict_num_decimal_digits( $totalG1Species,4);
			$totalG1Reads=&restrict_num_decimal_digits($totalG1Reads,4);
			
			
			#paired
			my %pairedG1guideStat=();
			$pairedG1guideStatRef=&g1Frac($ppgseqm);
			%pairedG1guideStat=%{$pairedG1guideStatRef};
			my $pairedG1Species=0;
			my $pairedG1Reads=0;
			my %pairedG1SpeciesC=();
			my %pairedG1ReadsC=();
			my %pairedG1SpeciesF=();
			my %pairedG1ReadsF=();
			foreach my $b (keys %pairedG1guideStat)
			{
				$pairedG1SpeciesC{$b}=scalar (keys %{$pairedG1guideStat{$b}});
				$pairedG1Species+=$pairedG1SpeciesC{$b};
				map {$pairedG1ReadsC{$b}+=$_} values %{$pairedG1guideStat{$b}};
				$pairedG1Reads+=$pairedG1ReadsC{$b};
			}
		
			foreach my $b(keys %pairedG1SpeciesC )
			{
				$pairedG1SpeciesF{$b}=$pairedG1SpeciesC{$b}/$pairedG1Species;
				$pairedG1SpeciesF{$b}=&restrict_num_decimal_digits($pairedG1SpeciesF{$b},4);
				$pairedG1ReadsF{$b}=$pairedG1ReadsC{$b}/$pairedG1Reads;
				$pairedG1ReadsF{$b}=&restrict_num_decimal_digits($pairedG1ReadsF{$b},4);
			}
			$pairedG1Species=&restrict_num_decimal_digits( $pairedG1Species,4);
			$pairedG1Reads=&restrict_num_decimal_digits($pairedG1Reads,4);
			
			
			print PPUAFRACTION "10\tg1totalByfile\tg1\t$pairedG1Species\t$totalG1Species\t$pairedG1Reads\t$totalG1Reads\n";
			foreach my $b(@bases)
			{
				print PPUAFRACTION "10\t$b\tg1Byfile\t$pairedG1SpeciesF{$b}\t$totalG1SpeciesF{$b}\t$pairedG1ReadsF{$b}\t$totalG1ReadsF{$b}\n";
			}
			
			$pptseq="$OUTDIR/$targetStrandFile.$basep.prefix.T.ppseq";
			$seqFile="$OUTDIR/$targetStrandFile.seq";
			$pptseqm="$OUTDIR/$targetStrandFile.$basep.prefix.T.ppseq.reads";
			`match.pl $pptseq $seqFile >$pptseqm`;
			
			my %totalT10guideStat=();
			
			$totalT10guideStatRef=&t10Frac($seqFile);
			%totalT10guideStat=%{$totalT10guideStatRef};
			my $totalT10Species=0;
			my $totalT10Reads=0;
			my %totalT10SpeciesC=();
			my %totalT10ReadsC=();
			my %totalT10SpeciesF=();
			my %totalT10ReadsF=();
			foreach my $b (keys %totalT10guideStat)
			{
				$totalT10SpeciesC{$b}=scalar (keys %{$totalT10guideStat{$b}});
				$totalT10Species+=$totalT10SpeciesC{$b};
				map {$totalT10ReadsC{$b}+=$_} values %{$totalT10guideStat{$b}};
				$totalT10Reads+=$totalT10ReadsC{$b};
			}
		
			foreach my $b(keys %totalT10SpeciesC )
			{
				$totalT10SpeciesF{$b}=$totalT10SpeciesC{$b}/$totalT10Species;
				$totalT10SpeciesF{$b}=&restrict_num_decimal_digits($totalT10SpeciesF{$b},4);
				$totalT10ReadsF{$b}=$totalT10ReadsC{$b}/$totalT10Reads;
				$totalT10ReadsF{$b}=&restrict_num_decimal_digits($totalT10ReadsF{$b},4);
			}
			$totalT10Species=&restrict_num_decimal_digits( $totalT10Species,4);
			$totalT10Reads=&restrict_num_decimal_digits($totalT10Reads,4);
			
			#paired
			my %pairedT10guideStat=();
			$pairedT10guideStatRef=&t10Frac($pptseqm);
			%pairedT10guideStat=%{$pairedT10guideStatRef};
			my $pairedT10Species=0;
			my $pairedT10Reads=0;
			my %pairedT10SpeciesC=();
			my %pairedT10ReadsC=();
			my %pairedT10SpeciesF=();
			my %pairedT10ReadsF=();
			foreach my $b (keys %pairedT10guideStat)
			{
				$pairedT10SpeciesC{$b}=scalar (keys %{$pairedT10guideStat{$b}});
				$pairedT10Species+=$pairedT10SpeciesC{$b};
				map {$pairedT10ReadsC{$b}+=$_} values %{$pairedT10guideStat{$b}};
				$pairedT10Reads+=$pairedT10ReadsC{$b};
			}
		
			foreach my $b(keys %pairedT10SpeciesC )
			{
				$pairedT10SpeciesF{$b}=$pairedT10SpeciesC{$b}/$pairedT10Species;
				$pairedT10SpeciesF{$b}=&restrict_num_decimal_digits($pairedT10SpeciesF{$b},4);
				$pairedT10ReadsF{$b}=$pairedT10ReadsC{$b}/$pairedT10Reads;
				$pairedT10ReadsF{$b}=&restrict_num_decimal_digits($pairedT10ReadsF{$b},4);
			}
			$pairedT10Species=&restrict_num_decimal_digits( $pairedT10Species,4);
			$pairedT10Reads=&restrict_num_decimal_digits($pairedT10Reads,4);
			
			
			print PPUAFRACTION "10\tT10totalByfile\tT10\t$pairedT10Species\t$totalT10Species\t$pairedT10Reads\t$totalT10Reads\n";
			foreach my $b(@bases)
			{
				print PPUAFRACTION "10\t$b\tT10Byfile\t$pairedT10SpeciesF{$b}\t$totalT10SpeciesF{$b}\t$pairedT10ReadsF{$b}\t$totalT10ReadsF{$b}\n";
			}
		}
	   
	   #total Ping-Pong
	   
		$score{$n}=0 if (!exists $score{$n});
		my $m=$n+1;
	    print PPSCORE "$m\t$score{$n}\n";
	    $count_N++ if ($score{$n}>0);
	   
	
	   
	   #Ping-Pong score according to different G1T10 pairs
	   #for matched pairs, cis only  
		foreach my $p (@matchedpairs)
		{
			$allPairReads{$p}{$n}=0 if (!exists $allPairReads{$p}{$n});
			my $n_of_allPairReads=$allPairReads{$p}{$n};
			$n_of_allPairReads=&restrict_num_decimal_digits($n_of_allPairReads,3);
			
			
			$n_of_species=scalar (keys %{$species{$p}{$n}});
			my $n_of_species_cor=0;
			map {$n_of_species_cor+=$_} values %{$species{$p}{$n}} ;
		     
			my $n_of_cisPairSpecies=0;
			$n_of_cisPairSpecies=scalar (keys %{$cisPairSpecies{$p}{$n}});					    
			my $n_of_cisPairSpecies_cor=0;
			map {$n_of_cisPairSpecies_cor+=$_} values %{$cisPairSpecies{$p}{$n}} ;					     
			my $n_of_cisPairReads=0;
			map {$n_of_cisPairReads+=$_} values %{$cisPairReads{$p}{$n}} ;
			$n_of_cisPairReads=&restrict_num_decimal_digits($n_of_cisPairReads,3);
			
			print PPSCOREUA "$m\tcis\t$p\t$n_of_cisPairSpecies\t$n_of_cisPairSpecies_cor\t$n_of_cisPairReads\t$n_of_species\t$n_of_species_cor\t$n_of_allPairReads\n";


	   }
	     #for all pairs, trans only 
		foreach my $p (@pairs)
		{
		     $allPairReads{$p}{$n}=0 if (!exists $allPairReads{$p}{$n});
		     my $n_of_allPairReads=$allPairReads{$p}{$n};
			 $n_of_allPairReads=&restrict_num_decimal_digits($n_of_allPairReads,3);
		     $n_of_species=scalar (keys %{$species{$p}{$n}});
		     my $n_of_species_cor=0;
			 map {$n_of_species_cor+=$_} values %{$species{$p}{$n}} ;
		     
		     my $n_of_transPairSpecies=0;					     					   
		     $n_of_transPairSpecies=scalar (keys %{$transPairSpecies{$p}{$n}});
		     my $n_of_transPairSpecies_cor=0;
		     map {$n_of_transPairSpecies_cor+=$_} values %{$transPairSpecies{$p}{$n}};					     
		     my $n_of_transPairReads=0;
		     map {$n_of_transPairReads+=$_} values %{$transPairReads{$p}{$n}} ;
		     $n_of_transPairReads=&restrict_num_decimal_digits($n_of_transPairReads,3);
		     
		    print PPSCOREUA "$m\ttrans\t$p\t$n_of_transPairSpecies\t$n_of_transPairSpecies_cor\t$n_of_transPairReads\t$n_of_species\t$n_of_species_cor\t$n_of_allPairReads\n";

		     $count_N0{$p}++ if ($n_of_transPairSpecies>0);
	     }
	     
	     
	     
	     
	}#n=1..20
	   $X=$score{9}; delete $score{9};
	   $std=&standard_deviation(values %score);
	   $m=&mean(values %score);
	   if ($std>0 && $count_N>=5) { $Z=($X-$m)/$std;} else {$Z=-10;}
	   $Z=&restrict_num_decimal_digits($Z,3);
	   $X=&restrict_num_decimal_digits($X,3);
	   $std=&restrict_num_decimal_digits($std,3);
	   $m=&restrict_num_decimal_digits($m,3);
	   print ZSCORE "$guideStrandFile\-$targetStrandFile\tall\tpp8\t$wsize\t$basep\t$Z\t$X\t$m\t$std\n";
	   #print "$guideStrandFile\-$targetStrandFile\t$Z\t";
	   
	   #Z-score for pp6
	  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%pp6cisPairSpecies,\%pp6cisPair10Species,\%pp6cisPairReads,$count_N);
	    #how to normalize $X0{$p}?
	    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tcisAll\tpp6\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
	   
	   #Z-score for individual cispairs; 
	   foreach my $p (@matchedpairs)
	   {
	  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%{$cisPairSpecies{$p}},\%{$cisPair10Species{$p}},\%{$cisPairReads{$p}},$count_N0{$p});
	    #how to normalize $X0{$p}?
	    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tcis\t$p\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
	   }
	   
	   #Z-score for individual transpairs; by species irrespective of coordinates	   
	   foreach my $p (@pairs)
	   {
	  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%{$transPairSpecies{$p}},\%{$transPair10Species{$p}},\%{$transPairReads{$p}},$count_N0{$p});
	    #how to normalize $X0{$p}?
	    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\ttrans\t$p\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
	   }
		
		#Z-score for all trans pairs
		my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%transallPairSpecies,\%transallPair10Species,\%transallPairReads,$count_N);
	    #how to normalize $X0{$p}?
	    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\ttransAll\tpp8\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n";	   

		
		#Z-score for pp8
	  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%pp8allPairSpecies,\%pp8allPair10Species,\%pp8allPairReads,$count_N);
	    #how to normalize $X0{$p}?
	    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tall\tpp8\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n";	   
		

	  
	close(PPTSEQ);  
	close(PPGSEQ);
	close(PPUAFRACTION);
	close(PPSCOREUA);
	close(PPSCORE);
	close(ZSCOREUA);
	close(ZSCORE);  	
}


sub ZscoreCal
{
		my ($transPairSpeciesRef,$transPair10SpeciesRef,$transPairReadsRef,$count)=@_;
		my ($Z0,$Z1,$Z2,$P0,$P1,$P2,$X0,$X1,$X2,$S0,$S1,$S2,$m0,$m1,$m2,$std0,$std1,$std2)=(0,0,0,0,0,0,0,0,0,0,0,0);
		
		#X0
		#$X0=scalar (keys %{$transPairSpeciesRef->{9}}); #by species irrespective of coordinates
	    #map {$X1+=$_} (values %{$transPairSpeciesRef->{9}}); #by species according to coordinates
	    #map {$X2+=$_} (values %{$transPairReadsRef->{9}});  #Z-score for all transpairs; by reads
	    
	    
	    
	    
	    
	    #for validation only
	    #$n_of_species=scalar (keys %{$transPair10Species}); #total number of species irrespective of coordinates
	    #my $n_of_species_cor=0;
	    #map {$n_of_species_cor+=$_} (values %{$transPair10Species});#total number of species respective of coordinates
	    #my $n_of_reads=0;
	    #map {$n_of_reads+=$_} (values %{$transPairReadsRef->{9}}); ###why it is not equal to $X2? pls calculate it before it's deleted!!
	    
	    #my %transPairSpecies9=%{$transPairSpeciesRef->{9}};
	    #my %transPairReads9=%{$transPairSpeciesRef->{9}};
	    #delete $transPairSpeciesRef->{9}; ##???
	    #delete $transPairReadsRef->{9};###???
	    
	    my @numOfSpecies=();
	    my @numOfSpeciesCor=();
	    my @numOfReads=();
	    
	    for(my $i=0;$i<$wsize;$i++)
	    {
	    	my $n1=scalar (keys %{$transPairSpeciesRef->{$i}});
	    	my $n2=scalar (keys %{$transPairReadsRef->{$i}});
	    	
	    	if($n1>0) #$transPairSpeciesRef->{9} was deleted
	    	{
	    		push @numOfSpecies, scalar (keys %{$transPairSpeciesRef->{$i}}) ;
	    		my $XnSpecies=0;
	    		map { $XnSpecies+=$_ } (values %{$transPairSpeciesRef->{$i}});
	    		push @numOfSpeciesCor, $XnSpecies ;
	    	}
			if($n2>0) #$transPairReadsRef->{9} was deleted
			{
				my $XnReads=0;
				map { $XnReads+=$_ } (values %{$transPairReadsRef->{$i}}); 
				push @numOfReads, $XnReads ;
			}				    	
	    }
	    $X0=$numOfSpecies[9];
	    $X1=$numOfSpeciesCor[9];
	    $X2=$numOfReads[9];
	    
	    $S0=sum(@numOfSpecies);
	    $S1=sum(@numOfSpeciesCor);
	    $S2=sum(@numOfReads);
	    
	    if($S0!=0)
	    {
	    	$P0=$X0/$S0;
	    }
	    else
	    {
	    	$P0=0;
	    	$X0=0;
	    }
	    
	    if($S1!=0)
	    {
	    	$P1=$X1/$S1;
	    }
	    else
	    {
	    	$P1=0;
	    	$X1=0;
	    }
	    
	    if($S2!=0)
	    {
	    	$P2=$X2/$S2;
	    }
	    else
	    {
	    	$P2=0;
	    	$X2=0;
	    }
	    
	 
	    
	    
	    
	    splice(@numOfSpecies, 9, 1);
	    splice(@numOfSpeciesCor, 9, 1);
	    splice(@numOfReads, 9, 1);
	    
	    if($S0!=0 && $S1!=0)
	    {
	    
		    $std0=&standard_deviation(@numOfSpecies);
		    $std1=&standard_deviation(@numOfSpeciesCor);
		    $std2=&standard_deviation(@numOfReads);
		    
		    
		    
		    #to prove that $transPairReadsRef->{9} was deleted successfully
		    #my $temp1=$#numOfSpecies+1;
		    #my $temp2=$#numOfSpeciesCor+1;
		    #my $temp3=scalar (@numOfReads);
		    
		    
		    $m0=&mean(@numOfSpecies);
		    $m1=&mean(@numOfSpeciesCor);
		    $m2=&mean(@numOfReads);
	    }
	    else
	    {
	    	$std0=0;
	    	$std1=0;
	    	$std2=0;
	    	$m0=0;
	    	$m1=0;
	    	$m2=0;
	    	
	    }
	    
	    if ($std0>0 && $count>=5) { $Z0=($X0-$m0)/$std0;} else {$Z0=-10;}#by species irrespective of coordinates
	    if ($std1>0 && $count>=5) { $Z1=($X1-$m1)/$std1;} else {$Z1=-10;}#by species according to coordinates
	    if ($std2>0 && $count>=5) { $Z2=($X2-$m2)/$std2;} else {$Z2=-10;}
	    
	    $Z0=&restrict_num_decimal_digits($Z0,3);
	    $Z1=&restrict_num_decimal_digits($Z1,3);
	    $Z2=&restrict_num_decimal_digits($Z2,3);
	    $X0=&restrict_num_decimal_digits($X0,3);
	    $X1=&restrict_num_decimal_digits($X1,3);
	    $X2=&restrict_num_decimal_digits($X2,3);
	    $m0=&restrict_num_decimal_digits($m0,3);
	    $m1=&restrict_num_decimal_digits($m1,3);
	    $m2=&restrict_num_decimal_digits($m2,3);
	    $std0=&restrict_num_decimal_digits($std0,3);
	    $std1=&restrict_num_decimal_digits($std1,3);
	    $std2=&restrict_num_decimal_digits($std2,3);
	    $P0=&restrict_num_decimal_digits($P0,4);
	    $P1=&restrict_num_decimal_digits($P1,4);
	    $P2=&restrict_num_decimal_digits($P2,4);
	    
	    
	    
	    
	    #%{$transPairSpeciesRef->{9}}=%transPairSpecies9;
	    #%{$transPairReadsRef->{9}}=%transPairReads9;
	    
	    return($Z0,$Z1,$Z2,$P0,$P1,$P2,$X0,$X1,$X2,$m0,$m1,$m2,$std0,$std1,$std2);
}

#remove intermediate files
#`rm *.ebwt`;
#`rm *.bowtie.out`;
sub g1Frac
{
	my $input= shift ;
	#input format
	#seq reads
	my %g1Stat=();
	open IN, $input;
	while(my $line=<IN>)
	{
		my($seq,$readn)=split(/\t/,$line);
		my $g1=substr($seq,0,1);
		$g1Stat{$g1}{$seq}+=$readn;

	}
	close(IN);
	return \%g1Stat;
}

sub t10Frac
{
	my $input= shift ;
	#input format
	#seq reads
	my %t10Stat=();
	open IN, $input;
	while(my $line=<IN>)
	{
		my($seq,$readn)=split(/\t/,$line);
		my $g1=substr($seq,9,1);
		$t10Stat{$g1}{$seq}+=$readn;

	}
	close(IN);
	return \%t10Stat;
}



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
        print "inputfile1 inputfile2 n type(fly|bombyx) outdir indexflag\n\t";
        print "-i  <inputfile1>\n\t";
		print "-j  <inputfile2[can be as same as inputfile1]>\n\t";
		print "-o  <outputdir>\n\t";
		print "-b  <bowtie index outputdir>\n\t";
		print "-n  <number of inputfiles>\n\t";
		print "-s  <species name[fly|bombyx]>\n\t";
		print "-d  <flag of index[0|1]>\n\t";
		print "-f  <flag of index[bed|normbed]>\n\t";
		print "-w  <background windowsize>\n\t";
		print "-p  <the length of prefix>\n\t";
		print "-a  <fasta file of the genome>\n\t";
		print "-m  <mapping output dir>\n\t";
		print "-q  <query seq output dir>\n\t";
        print "This perl script is count the frequency of 10A irrespective of 1U\n";
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
                elsif($next_arg eq "-n"){ $parameters->{num} = shift(@ARGV); }
                elsif($next_arg eq "-s"){ $parameters->{species} = shift(@ARGV); }
                elsif($next_arg eq "-o"){ $parameters->{outdir} = shift(@ARGV); }
                elsif($next_arg eq "-b"){ $parameters->{indexoutdir} = shift(@ARGV); }
                elsif($next_arg eq "-d"){ $parameters->{indexflag} = shift(@ARGV); }
                elsif($next_arg eq "-f"){ $parameters->{format} = shift(@ARGV); }
                elsif($next_arg eq "-w"){ $parameters->{winsize}= shift(@ARGV); }
				elsif($next_arg eq "-p"){ $parameters->{complementarity}= shift(@ARGV); }
				elsif($next_arg eq "-a"){ $parameters->{fa} = shift(@ARGV); }
				elsif($next_arg eq "-m"){ $parameters->{mappingoutdir}= shift(@ARGV); }
				elsif($next_arg eq "-q"){ $parameters->{queryseqoutdir} = shift(@ARGV); }

                else{ print "Invalid argument: $next_arg"; usage(); }
        }
}
