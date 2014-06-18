#!/usr/bin/perl

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "restrict_digts.pm";
require "Jia.pm";
use File::Basename;
use Compress::Zlib;
#use List::Util qw(sum);
use PDL;  #it's needed for initialize, it has sum functio inside, clash with the one in Util
use PDL::Char;
$PDL::SHARE = $PDL::SHARE; 

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

#03/05
#save memory
#remove %allPairReads, %species,%species10
#report each piRNA species which have Ping-Pong pairs

#03/08
#For guide or target that are multi-mappers, if one out of all the possible combinations of mapping coordinates is cis, this pair is considered cis.

#03/09
#add mapping output dir and query seq output dir

#03/10
#store each piRNA species when they are collapsed into the same $basep prefix 
#calculate the number of prefix species, the number of species for pp10

#03/13/2014
#don't if cis-pairs are assigned correctly, try to save some trans-pairs from last version
#do I need to reconsider NTM?

#5/15/2014
#report Ping-Pong pairs

#06/16/2014
#change the cis pairs to pure cis pairs(pp6);
#new method to check how much complementarity is required
#ask for full complementarity at position 2-10
#check if the rest are matched(1) or not matched(0), species, or reads(scale by normalized reads)

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

 my $supLen=23-$basep;



open LOG, ">${OUTDIR}/LOG.txt";


my @matchedpairs=("AT","TA","GC","CG");
my @unmatchedpairs=("AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");
my @pairs=("AT","TA","GC","CG","AA","AC","AG","CA","CC","CT","GA","GG","GT","TC","TG","TT");



#my %targetpf=();
#my %guidepf=(); 

my %targetpfsplit=();
my %guidepfsplit=();

my %total=();
my %genome=();
#store genome
if($spe eq "fly")
{
	open IN, $fastafile;
	while(my $line=<IN>)
	{
	   if ($line=~/>(.+) type/)
	   {
	      $chr="chr$1";
	   }
	   else
	   {
	      chomp $line;
	      $genome{$chr}=$line;
	   }
	}
	close(IN);
}
elsif($spe eq "bombyx")
{

	open IN, $fastafile or die "Fail to open $fastafile: $!";
	while(my $line=<IN>)
	{
	   if ($line=~/>(.+)\s*\//) #this is specific for the case: >nscaf100 /length=4083 /lengthwogaps=4073
	   {
	      $chr="$1";
	      @c=split(/ /,$chr);
	   }
	   else
	   {
	      chomp $line;
	      $genome{$c[0]}=$line;
	   }
	}
	close(IN);
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
#memory


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
	    
	    $file=$name;
	    
	    &InputFileProcessing($inputfiles[$i],$file);
	    

		
		
		#if ($total{$file}>10)
		#{
			#query seq file 
			$seqFile="$QOUTDIR/$file.$basep.seq";
			if( ! -s $seqFile )#test the existence of file
			{
				open OUT, ">$seqFile";
				foreach my $prefix (keys %{$guidepf{$file}})
				{
					print OUT "$prefix\t$guidepf{$file}{$prefix}\n" if (length($prefix)==$basep);
				}
				close(OUT);
			}
			#bowtie index file
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
					close(OUT);
					`bowtie-build $fa $indexb && rm $fa`;
				}
	
			}#for loop of the n
		#}#if the total
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
		#my $memnow=qx{ `grep -i VmSize /proc/$$/status` };
		#print LOG "the memory used for preprocessing is: $memnow";
   


		#if ($total{$file1}<10 || $total{$file2}<10) {print PPZ "$file2-$file1\t-10\n";} #only when consider per cluster or per transposon family
		#else
		#{
			&PingPongProcessing($file2,$file1);

			if($file1 ne $file2 ) #added on 11/14/2013
			{   				
			   &PingPongProcessing($file1,$file2); 

			}#when file1 and file2 are the same the second iteration is skipped
		#}#ifelse: total reads>10 

	}#j
}#i
close(PPZ);
close(LOG);



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
		my $len;
		if($fileFormat eq "normbed")
		{ 
			($chr,$bedstart,$bedend,$strand,$seq,$reads,$ntm)= split(/\t/,$line);
			$len=$bedend-$bedstart+1;
		}
		if($fileFormat eq "bed")
		{
			($chr,$bedstart,$bedend,$seq,$ntmreads,$strand)= split(/\t/,$line);
			$reads=$ntmreads;
			$ntm=1;
			$len=$bedend-$bedstart;
		}
		
		if($inputfile=~/SRA/){next if (length($seq)>29 || length($seq)<23);}
		next if (/data/);
		
		$total{$file}+=$reads/$ntm; #no nta in norm.bed format
		#store chr, 5'end and strand information separately for each query 20nt prefix, the populations to find the guide
		my $fiveend=0;
		my $dnaseq="";
		my $supplementseq="";
		my $piRNA="";	      
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
			#the piRNA species
			$piRNA=substr($genome{$chr},$fiveend,$len); #this is from genome strand (+ strand)
			#the prefix of piRNA species
			$dnaseq=substr($piRNA,0,$basep);
			
			#the supplemental region of piRNAs
			#my $suppLen=23-$basep;
			$supplementseq=substr($piRNA,$basep,$supLen);
					
			#store the seq of guide 20nt prefix only; for faster extract the reads number later
			$guidepf{$file}{$dnaseq}+=$reads/$ntm;
			
			#the coordinates of each piRNA species
			$guidepfsplit{$file}{$dnaseq}{"$piRNA,$supplementseq"}{"$chr,$fiveend,$strand"}+=$reads/$ntm; #become 0-based from norm.bed format


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
			my $seqstart=$fiveend-$len; #fix this bug on 03-14-2014; this is the full length of piRNA species, not the same as prefix
     		my $seqtemp=substr($genome{$chr},$seqstart,$len); #this is the genomic strand (+ strand)
     		$piRNA=&revfa($seqtemp);
			$dnaseq=substr($piRNA,0,$basep);
			#my $suppLen=23-$basep;
			$supplementseq=substr($piRNA,$basep,$supLen);			
			#store the seq of guide 20nt prefix only; for faster extract the reads number later
			$guidepf{$file}{$dnaseq}+=$reads/$ntm;
			$guidepfsplit{$file}{$dnaseq}{"$piRNA,$supplementseq"}{"$chr,$fiveend,$strand"}+=$reads/$ntm; #become 0-based from norm.bed format

	  	}

      	for (my $n=0;$n<$wsize;$n++)
      	{
      		my $start=0;
      		my $supStart=0;
      		my $fiveend=0;
        	if ($strand eq '+') #target strand information
         	{
         		
         		if($fileFormat eq "bed")
         		{
	            	$start=$bedstart+$n+1-$basep; # $bedstart is 0 based and $start is 0 based, the intermediate end $bedstart+$n+1 is open
	            	$supStart=$bedstart+$n+1-23;
	            	$fiveend=$start+$basep; #open
         		}
         		if($fileFormat eq "normbed")
         		{
         			$start=$bedstart+$n-$basep;
         			$supStart=$bedstart+$n-23;
         			$fiveend=$start+$basep;#convert to bed, open
         		}
	            
	            
	            
	            #supplemental region
	           
	            my $supStr=substr($genome{$chr},$supStart,$supLen);
	            
	            if(length($supStr) == $supLen) #exclude boundary reads
	            {
	            	
	            my $str=substr($genome{$chr},$start,$basep); #substr function is 0 based
	            $str=&revfa($str);
	            $targetpf{$file}{$n}{$str}+=$reads/$ntm;
	            	
	            $subStr=&revfa($supStr);
	            
	            #store chr, 5'end and strand information separately for each guide 16nt prefix
	            my $tstrand="-";
	            $targetpfsplit{$file}{$n}{$str}{"$piRNA,$subStr"}{"$chr,$fiveend,$tstrand"}+=$reads/$ntm; #store the strand information for guide strand
	            #my $indexStart=$start;
	            }
	            	            
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
	            
	            $supStart=$start+$basep;
	            my $supStr=substr($genome{$chr},$supStart,$supLen);
	            if(length($supStr) == $supLen) #exclude boundary reads
	            {
	            	my $str=substr($genome{$chr},$start,$basep);
	            	$targetpf{$file}{$n}{$str}+=$reads/$ntm;
	            	my $tstrand="+";
	            	#store chr, 5'end and strand information separately for each guide 16nt prefix	
	            	$targetpfsplit{$file}{$n}{$str}{"$piRNA,$subStr"}{"$chr,$fiveend,$tstrand"}+=$reads/$ntm;
	            }
	            	            
        	}#ifelse
      	}#for
	} #while
	$gz->gzclose();
}


sub PingPongProcessing
{

	my ($guideStrandFile,	$targetStrandFile)=@_;		


	my $count_N=0;
	my %count_N0=();

	
	my %cisPairSpecies=();
	my %cisPairReads=();
	
	my %ambiguousPairSpecies=();	
	my %ambiguousPairReads=();
	
	my %transPairSpecies=();
	my %transPairReads=();
	

	
	my %transPairSuppSpeciesTotal=();
	my %transPairSuppReadsTotal=();	

	open ZSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.zscore.out";
	open PPSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.pp";
	open PPSEQPAIR, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.ppseq.txt";
	open PPSEQSUPPVECTOR, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.ppsub.vector.scale.txt";
	
	#open PPUAFRACTION, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.base.fraction.txt";
	print PPSEQPAIR "pairmode\tguidepiRNAs\tguidepiRNAsReads\ttargetpiRNAs\ttargetpiRNAsReads\n";

	
	foreach ($n=0;$n<$wsize;$n++)
	{
		
		my %cisPairSuppSpecies=();
		my %cisPairSuppReads=();
		
		my %ambiguousPairSuppSpecies=();
		my %ambiguousPairSuppReads=();
		
		my %transPairSuppSpecies=();
		my %transPairSuppReads=();

		# file1 as ref
		$indexb="$BOUTDIR/$targetStrandFile.$basep.$n";
		$seqFile="$QOUTDIR/$guideStrandFile.$basep.seq";
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
		   	$NTM{$l[2]}++; #need to think about why only query sequences are nomalized to NTM, but not indexes
	   	}
	   	close(IN);


	   	open IN, "$MOUTDIR/$guideStrandFile.$targetStrandFile.$basep.$n.bowtie.out";
	   	while(my $line=<IN>)
	   	{
	      	chomp $line;
	      	@l=split(/\t/,$line);
	
	      	my $nGcor=scalar (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}}); #total piRNA species for this prefix from guide strand
			
			my $nGcorTotalReads=0; #total reads for this prefix from guide strand
			foreach my $piQuery (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}})
			{
				my $nGcorReads=0;
				map {$nGcorReads+=$_} values %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}};
				$nGcorTotalReads+=$nGcorReads;
			}
			
		    my $nTcor=scalar (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}}); #total piRNA species for this prefix from target strand		
		    
			my $nTcorTotalReads=0; #total reads for this prefix from target strand
			foreach my $piTargetIndex (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}} )
			{
				my $nTcorReads=0;
				map {$nTcorReads+=$_} values %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}};
				$nTcorTotalReads+=$nTcorReads;
			}
		    
		    my $nnGcorTcor=$nGcor*$nTcor; #total species pairs for this prefix
	      	my $gttotal=$nGcorTotalReads*$nTcorTotalReads; #total read pairs for this prefix
	      	if ($l[3] eq "")
	      	{
	      	
	      	   $g_0_nt=substr($l[2],0,1); $t_9_nt=&revfa($g_0_nt);  ##here are different from pp8_q2_ww1.pl
		       #targetpf index; guidepf seq
			   #how many of species start with U?
  		        		       		     		      		       
	       	   	#my %cisFlag=(); #this cisFlag is to mark if there is a cis pair among all the possible pairs(by coordinates) among a paired species	    		       
				foreach my $piQuery (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}} )
			    {
			    	
			    	my ($piRNAGuide, $piGuideSuppSeq)=split(/,/,$piQuery);
					my $guideQuerySpecies=0;
					$guideQuerySpecies=scalar (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}}); #for this guide piRNA, the number of mapping locus(genomic coordinates) 
					my $guideQueryReads=0;
					map {$guideQueryReads+=$_} values %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}};#for this guide piRNA, the number of reads for all mapping locus 
					#guide is the query, target is the index
      				my $guideQueryReadsNorm=$guideQueryReads/$NTM{$l[2]};	
      				
					foreach my $piTargetIndex (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}} ) #indexed piRNA prefix; iterate each piRNA
					{	

						#my $targetpiGuideSpecies=0;	
						#$targetpiGuideSpecies=scalar (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}});
						
						my ($targetpiRNA,$targetSuppSeq)=split(/,/,$piTargetIndex);
						my $ref=$targetSuppSeq;
	      				my $source = PDL::Char->new($ref);
	      				my $match = PDL::Char->new($piGuideSuppSeq);     							
						my $diff = $match == $source;
						my $diffstr=join('',$diff->list);
						
						my $targetpiGuideReads=0;
						map {$targetpiGuideReads+=$_} values %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}}; # the total reads for this piRNA
						
						my $cisRecordFlag=0; #flag of cispair
						
						foreach my $record (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}}) #each mapping locus for guide piRNA
						{
				       		my ($chr,$gfiveend,$gstrand)=split(/,/,$record);
				       		
				       		my $guideQueryCorReads=$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}{$record};
				       		my $guideQueryCorReadsNorm=$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}{$record}/$NTM{$l[2]};
				       		
				       		###note: how to define cistargets more accurately?
				       		###in addition to excat match, what if just several nucleotides away?

				       		my $tfiveend=$gfiveend;
				       		my $tstrand=$gstrand;
							my $targetpiGuideCorReads=$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}{"$chr,$tfiveend,$tstrand"} ;
							# for each mapping locus for this guide piRNA, if there is a cis pair from target piRNA coordinate pool	
					       	if($targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}{"$chr,$tfiveend,$tstrand"}) #here it checks all coordinates associated with $piTargetIndex
					       	{
					       		my $nGpiCor=scalar (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}});
      							if( $nGpiCor>1)
					       		{
						       		$cisPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=1/$NTM{$l[2]}/2; #cis pair species(in terms of piRNA species) must only have one by coordinate definition; but for different record, it has different cis pair
									$cisPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guideQueryCorReads*$targetpiGuideCorReads/$NTM{$l[2]};		       			      							
	      							##print out paired piRNA species      							
	      							print PPSEQPAIR "cis\t$piQuery\t$guideQueryCorReadsNorm\t$piTargetIndex\t$targetpiGuideCorReads\t$diffstr\n" if ($n==9);
	      								      								     							
	     							$cisPairSuppSpecies{$g_0_nt.$t_9_nt}{$diffstr}+=1/$NTM{$l[2]}/2;
	      							$cisPairSuppReads{$g_0_nt.$t_9_nt}{$diffstr}+=$guideQueryCorReads*$targetpiGuideCorReads/$NTM{$l[2]};
	      								      							
	      							#ambiguous pairs if there are other mapping loci
      								$ambiguousPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=1/$NTM{$l[2]}/2;
      								$ambiguousPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=($guideQueryReads-$guideQueryCorReads)*($targetpiGuideReads-$targetpiGuideCorReads)/$NTM{$l[2]};
      								
      								my $ambiguousGuideQueryCorReadsNorm=$guideQueryReadsNorm-$guideQueryCorReadsNorm;
      								my $ambiguousTargetpiGuideCorReads=$targetpiGuideReads-$targetpiGuideCorReads;
      								print PPSEQPAIR "ambiguous\t$piQuery\t$ambiguousGuideQueryCorReadsNorm\t$piTargetIndex\t$ambiguousTargetpiGuideCorReads\t$diffstr\n" if ($n==9);
      								
      								$ambiguousPairSuppSpecies{$g_0_nt.$t_9_nt}{$diffstr}+=1/$NTM{$l[2]}/2;
      								$ambiguousPairSuppReads{$g_0_nt.$t_9_nt}{$diffstr}+=($guideQueryReads-$guideQueryCorReads)*($targetpiGuideReads-$targetpiGuideCorReads)/$NTM{$l[2]};
      								
      							}
      							else
      							{
      								$cisPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=1/$NTM{$l[2]};
      								$cisPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guideQueryCorReads*$targetpiGuideCorReads/$NTM{$l[2]};		       			      							
	      							##print out paired piRNA species      							
	      							print PPSEQPAIR "cis\t$piQuery\t$guideQueryCorReadsNorm\t$piTargetIndex\t$targetpiGuideCorReads\t$diffstr\n" if ($n==9);
	      								     							
	     							$cisPairSuppSpecies{$g_0_nt.$t_9_nt}{$diffstr}+=1/$NTM{$l[2]}/2;
	      							$cisPairSuppReads{$g_0_nt.$t_9_nt}{$diffstr}+=$guideQueryCorReads*$targetpiGuideCorReads/$NTM{$l[2]};	
      							}
      							      						
								#$cisFlag{$piTargetIndex}+=1;
								$cisRecordFlag=1;
								last; #no need to check other coordinates (record) for this pair of piRNA species      					       					       			
					       	}#if cis coordinate exists
						} #record
						#for each query(a piRNA species), if there is a cismatch, then the whole species (no matter how many other possible pairs it has) is viewed as cispair
						if (! $cisRecordFlag) #after check all the records associated with piQuery
			       		{
			       			#trans PingPong pair in species
			       			$transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=1/$NTM{$l[2]};
			       			#trans PingPong pair in reads
			       			$transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guideQueryReads*$targetpiGuideReads/$NTM{$l[2]};
			       			print PPSEQPAIR "trans\t$piQuery\t$guideQueryReadsNorm\t$piTargetIndex\t$targetpiGuideReads\t$diffstr\n" if ($n==9);
			       			
			       			$transPairSuppSpecies{$g_0_nt.$t_9_nt}{$diffstr}+=1/$NTM{$l[2]};
			       			$transPairSuppReads{$g_0_nt.$t_9_nt}{$diffstr}+=$guideQueryReads*$targetpiGuideReads/$NTM{$l[2]};
			       			
			       			$transPairSuppSpeciesTotal{$diffstr}+=1/$NTM{$l[2]} if ($n==9);
			       			$transPairSuppReadsTotal{$diffstr}+=$guideQueryReads*$targetpiGuideReads/$NTM{$l[2]} if ($n==9);
			       			
			       		}
						
						#if(! $cisFlag{$piTargetIndex}) #so that no need to check for already cispaired piRNA species
						#{
						#}
					}#piGuideSpe
			    }#piSpecies
	      	}#perfect pair

	      	elsif ($l[3]=~/(\d+):(\w)>(\w)/)
	      	{

		       next if ($1!=0);  # allow 1mm at the 10th position of target strand
		       
		       $g_0_nt=$3;
		       $t_9_nt=&revfa($2);
		       	
		       #trans PingPong pair in species
		       $transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$nnGcorTcor/$NTM{$l[2]};
		       #trans PingPong pair in reads
		       $transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$gttotal/$NTM{$l[2]};
		       		       		       		       		      		   
		       	foreach my $piQuery (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}} )
			    {
					my ($piRNAGuide, $piGuideSuppSeq)=split(/,/,$piQuery);
					
					my $guideQueryReads=0;
					map {$guideQueryReads+=$_} values %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}};
					my $guideQueryReadsNorm=$guideQueryReads/$NTM{$l[2]};	
					foreach my $piTargetIndex (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}} )
					{	
						my ($targetpiRNA,$targetSuppSeq)=split(/,/,$piTargetIndex);
						my $ref=$targetSuppSeq;
	      				my $source = PDL::Char->new($ref);
	      				my $match = PDL::Char->new($piGuideSuppSeq);     							
						my $diff = $match == $source;
						my $diffstr=join('',$diff->list);
																	
						#my $targetpiGuideSpecies=0;	
						#$targetpiGuideSpecies=scalar (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}});
						my $targetpiGuideReads=0;
						map {$targetpiGuideReads+=$_} values %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}};
		       			print PPSEQPAIR "trans\t$piQuery\t$guideQueryReadsNorm\t$piTargetIndex\t$targetpiGuideReads\t$diffstr\n" if ($n==9);
		       			
		       			
		       			$transPairSuppSpecies{$g_0_nt.$t_9_nt}{$diffstr}+=1/$NTM{$l[2]};
			   			$transPairSuppReads{$g_0_nt.$t_9_nt}{$diffstr}+=$guideQueryReads*$targetpiGuideReads/$NTM{$l[2]};
			   			
			   			$transPairSuppSpeciesTotal{$diffstr}+=1/$NTM{$l[2]} if ($n==9);
			       		$transPairSuppReadsTotal{$diffstr}+=$guideQueryReads*$targetpiGuideReads/$NTM{$l[2]} if ($n==9);
			   			
					}
			    }
			}
		} #while
		close(IN);
        
		$m=$n+1;
	
			   #Ping-Pong score according to different G1T10 pairs
			   #for matched pairs, cis only  
				foreach my $p (@matchedpairs)
				{  
					my $n_of_cisPairSpecies=0;
					$n_of_cisPairSpecies=scalar (keys %{$cisPairSpecies{$p}{$n}});
					$n_of_cisPairSpecies=&restrict_num_decimal_digits($n_of_cisPairSpecies,3);
									    
					my $n_of_cisPairSpecies_cor=0;
					map {$n_of_cisPairSpecies_cor+=$_} values %{$cisPairSpecies{$p}{$n}} ;
					$n_of_cisPairSpecies_cor=&restrict_num_decimal_digits($n_of_cisPairSpecies_cor,3);
											     
					my $n_of_cisPairReads=0;
					map {$n_of_cisPairReads+=$_} values %{$cisPairReads{$p}{$n}} ;
					$n_of_cisPairReads=&restrict_num_decimal_digits($n_of_cisPairReads,3);

					print PPSCOREUA "$m\tcis\t$p\t$n_of_cisPairSpecies\t$n_of_cisPairSpecies_cor\t$n_of_cisPairReads\n";
					
					my $n_of_ambiguousPairSpecies=0;
					$n_of_ambiguousPairSpecies=scalar (keys %{$ambiguousPairSpecies{$p}{$n}});
					$n_of_ambiguousPairSpecies= &restrict_num_decimal_digits($n_of_ambiguousPairSpecies,3);
										    
					my $n_of_ambiguousPairSpecies_cor=0;
					map {$n_of_ambiguousPairSpecies_cor+=$_} values %{$ambiguousPairSpecies{$p}{$n}} ;
					$n_of_ambiguousPairSpecies_cor=&restrict_num_decimal_digits($n_of_ambiguousPairSpecies_cor,3);
											     
					my $n_of_ambiguousPairReads=0;
					map {$n_of_ambiguousPairReads+=$_} values %{$ambiguousPairReads{$p}{$n}} ;
					$n_of_ambiguousPairReads=&restrict_num_decimal_digits($n_of_ambiguousPairReads,3);

					print PPSCOREUA "$m\tambiguous\t$p\t$n_of_ambiguousPairSpecies\t$n_of_ambiguousPairSpecies_cor\t$n_of_ambiguousPairReads\n";
	
			   }
			   #Ping-Pong score according to different G1T10 pairs
			     #for all pairs, trans only 
				foreach my $p (@pairs)
				{
				     my $n_of_transPairSpecies=0;					     					   
				     $n_of_transPairSpecies=scalar (keys %{$transPairSpecies{$p}{$n}});
				     $n_of_transPairSpecies=&restrict_num_decimal_digits($n_of_transPairSpecies,3);
				     
				     my $n_of_transPairSpecies_cor=0;
				     map {$n_of_transPairSpecies_cor+=$_} values %{$transPairSpecies{$p}{$n}};
					 $n_of_transPairSpecies_cor=&restrict_num_decimal_digits($n_of_transPairSpecies_cor,3);

				     my $n_of_transPairReads=0;
				     map {$n_of_transPairReads+=$_} values %{$transPairReads{$p}{$n}} ;
				     $n_of_transPairReads=&restrict_num_decimal_digits($n_of_transPairReads,3);

				     print PPSCOREUA "$m\ttrans\t$p\t$n_of_transPairSpecies\t$n_of_transPairSpecies_cor\t$n_of_transPairReads\n";

				    $count_N0{$p}++ if ($n_of_transPairSpecies>0);
			     }

				 #supplemental pairing score according to different g1t10 pairs, no $n was recorded in hash
				 foreach my $p (@pairs)
				 {
				 	my %transPairIndiSuppSpecies= %{$transPairSuppSpecies{$p}};
				 	my %transPairIndiSuppReads= %{$transPairSuppReads{$p}};
				 	
				 	my ($scaledSpeciesRef,$scaledReadsRef)=&SuppComBitSum(\%transPairIndiSuppSpecies,\%transPairIndiSuppReads);
					for(my $position=0; $position< @{$scaledSpeciesRef};$position++)
					{
						my $supPos=$position+$basep+1;				
						print PPSEQSUPPVECTOR "transPair\t$p\t$basep\t$m\t$supPos\t$scaledSpeciesRef->[$position]\t$scaledReadsRef->[$position]\n";
					}
				 	
				 }
				foreach my $p (@matchedpairs)
				{
					my %cisPairIndiSuppSpecies= %{$cisPairSuppSpecies{$p}};
				 	my %cisPairIndiSuppReads= %{$cisPairSuppReads{$p}};
				 	
				 	my ($scaledSpeciesRef,$scaledReadsRef)=&SuppComBitSum(\%cisPairIndiSuppSpecies,\%cisPairIndiSuppReads);#
					for(my $position=0; $position< @{$scaledSpeciesRef};$position++)
					{
						my $supPos=$position+$basep+1;				
						print PPSEQSUPPVECTOR "cisPair\t$p\t$basep\t$m\t$supPos\t$scaledSpeciesRef->[$position]\t$scaledReadsRef->[$position]\n";
					}
					
					
					my %ambiguousPairIndiSuppSpecies= %{$ambiguousPairSuppSpecies{$p}};
				 	my %ambiguousPairIndiSuppReads= %{$ambiguousPairSuppReads{$p}};
				 	
				 	my ($scaledSpeciesRef,$scaledReadsRef)=&SuppComBitSum(\%ambiguousPairIndiSuppSpecies,\%ambiguousPairIndiSuppReads);#for n=10
					for(my $position=0; $position< @{$scaledSpeciesRef};$position++)
					{
						my $supPos=$position+$basep+1;				
						print PPSEQSUPPVECTOR "ambiguousPair\t$p\t$basep\t$m\t$supPos\t$scaledSpeciesRef->[$position]\t$scaledReadsRef->[$position]\n";
					}

				}
				 
				 


			}#n=1..16
	

			   	my %pp8allPairSpecies=();
				my %pp8allPairReads=();

				my %pp6cisPairSpecies=();
				my %pp6cisPairReads=();

				foreach my $p (keys %cisPairSpecies)
				{
					foreach my $num (keys %{$cisPairSpecies{$p}})
				    {
				       	foreach my $guide (keys %{$cisPairSpecies{$p}{$num}})
				       	{
				       		$pp6cisPairSpecies{$num}{$guide}+=$cisPairSpecies{$p}{$num}{$guide};
				       		$pp8allPairSpecies{$num}{$guide}+=$cisPairSpecies{$p}{$num}{$guide};
				       	}
				    }
				}
				foreach my $p (keys %cisPairReads)
				{
					foreach my $num (keys %{$cisPairReads{$p}})
				    {
				       	foreach my $guide (keys %{$cisPairReads{$p}{$num}})
				       	{
				       		$pp6cisPairReads{$num}{$guide}+=$cisPairReads{$p}{$num}{$guide};
				       		$pp8allPairReads{$num}{$guide}+=$cisPairReads{$p}{$num}{$guide};
				       	}
				    }
				}

			

				
				print ZSCOREUA "guide-target\tpairmode\tstatmode\twindowSize\tbasePairingext\tZscoreofprefixSpecies\tZscoreofSpecies\tZscoreofpairsofReads\tPercentageofprefixSpecies\tPercentageofSpecies\tPercentageofparisofReads\tnumofprefixSpeciesofpp10\tnumofSpeciesofpp10\tpairsofReadsofpp10\tmeanofprefixSpecies\tmeanofSpecies\tmeanofpairsofReads\tstdofprefixSpecies\tstdofSpecies\tstdofpairsofReads\n";		   	   

			   #Z-score for pp6
			  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%pp6cisPairSpecies,\%pp6cisPairReads);
			    #how to normalize $X0{$p}?
			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tcisAll\tcisAll4\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
			    undef %pp6cisPairSpecies;
			    undef %pp6cisPairReads;

			   #Z-score for individual cispairs; 
			   foreach my $p (@matchedpairs)
			   {
			  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%{$cisPairSpecies{$p}},\%{$cisPairReads{$p}});
			    #how to normalize $X0{$p}?
			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tcis\t$p\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
			   }

			   #Z-score for individual transpairs; by species, reads and by species irrespective of coordinates	   
			   foreach my $p (@pairs)
			   {
			  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%{$transPairSpecies{$p}},\%{$transPairReads{$p}});
			    #how to normalize $X0{$p}?
			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\ttrans\t$p\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
			   }
			   
			   #Z-score for individual ambiguouspairs; by species, reads and by species irrespective of coordinates	   
			   foreach my $p (@matchedpairs)
			   {
			  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%{$ambiguousPairSpecies{$p}},\%{$ambiguousPairReads{$p}});
			    #how to normalize $X0{$p}?
			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tambiguous\t$p\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
			   }


				my %transallPairSpecies=();
				my %transallPairReads=();

				foreach my $p (keys %transPairSpecies)
				{
					foreach my $num (keys %{$transPairSpecies{$p}})
				    {
				       	foreach my $guide (keys %{$transPairSpecies{$p}{$num}})
				       	{
				       		$transallPairSpecies{$num}{$guide}+=$transPairSpecies{$p}{$num}{$guide};
				       		$pp8allPairSpecies{$num}{$guide}+=$transPairSpecies{$p}{$num}{$guide};
				       	}
				    }
				}
				foreach my $p (keys %transPairReads)
				{
					foreach my $num (keys %{$transPairReads{$p}})
				    {
				       	foreach my $guide (keys %{$transPairReads{$p}{$num}})
				       	{
				       		$transallPairReads{$num}{$guide}+=$transPairReads{$p}{$num}{$guide};
				       		$pp8allPairReads{$num}{$guide}+=$transPairReads{$p}{$num}{$guide};
				       	}
				    }
				}
				
				#my $memnow=qx{ `grep -i VmSize /proc/$$/status` };
				#print LOG "the memory used for constructing transall and pp8 hashes during $guideStrandFile.$targetStrandFile processing is: $memnow";

				#Z-score for all trans pairs;
				my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%transallPairSpecies,\%transallPairReads);
			    #how to normalize $X0{$p}?
			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\ttransAll\ttransAll16\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n";	   
				undef %transallPairSpecies;
				undef %transallPairReads;

				#Z-score for pp8;
			  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%pp8allPairSpecies,\%pp8allPairReads);
			    #how to normalize $X0{$p}?
			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tall\tpp8\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n";	   

				undef  %pp8allPairSpecies;
				undef %pp8allPairReads;
				
				#score for supplemental region pairing if view seed region as 2-10
				my ($scaledSpeciesRef,$scaledReadsRef)=&SuppComBitSum(\%transPairSuppSpeciesTotal,\%transPairSuppReadsTotal);#for n=10
				for(my $position=0; $position< @{$scaledSpeciesRef};$position++)
				{
					my $supPos=$position+$basep+1;				
					print PPSEQSUPPVECTOR "transPair\tall\t$basep\t10\t$supPos\t$scaledSpeciesRef->[$position]\t$scaledReadsRef->[$position]\n";
				}
				#for(my $position=0; $position< @{$scaledReadsRef};$position++)
				#{
				#	my $supPos=$position+$basep+1;				
				#	print PPSEQSUPPVECTOR "reads\t$supPos\t$scaledReadsRef->[$position]\n";
				#}


			#close(PPUAFRACTION);
			close(PPSEQSUPPVECTOR);
			close(PPSEQPAIR);
			close(PPSCOREUA);
			close(ZSCOREUA);
 	
		}

		sub ZscoreCal
		{
				my ($transPairSpeciesRef,$transPairReadsRef)=@_;
				my ($Z0,$Z1,$Z2,$P0,$P1,$P2,$X0,$X1,$X2,$S0,$S1,$S2,$m0,$m1,$m2,$std0,$std1,$std2)=(0,0,0,0,0,0,0,0,0,0,0,0);

			    my @numOfSpecies=();
			    my @numOfSpeciesCor=();
			    my @numOfReads=();

			    for(my $n=0;$n<$wsize;$n++)
			    {
			    	my $n1=scalar (keys %{$transPairSpeciesRef->{$n}});
			    	my $n2=scalar (keys %{$transPairReadsRef->{$n}});

			    	if($n1>0) #$transPairSpeciesRef->{9} was deleted
			    	{
			    		push @numOfSpecies, scalar (keys %{$transPairSpeciesRef->{$n}}) ;
			    		my $XnSpecies=0;
			    		map { $XnSpecies+=$_ } (values %{$transPairSpeciesRef->{$n}});
			    		push @numOfSpeciesCor, $XnSpecies ;
			    	}
					else
					{
						push @numOfSpecies,0;
						push @numOfSpeciesCor,0;
					}
					if($n2>0) #$transPairReadsRef->{9} was deleted
					{
						my $XnReads=0;
						map { $XnReads+=$_ } (values %{$transPairReadsRef->{$n}}); 
						push @numOfReads, $XnReads ;
					}
					else
					{
						push @numOfReads,0;
					}				    	
			    }
			    $X0=$numOfSpecies[9];
			    $X1=$numOfSpeciesCor[9];
			    $X2=$numOfReads[9];

			    $S0=&arrSum(@numOfSpecies);
			    $S1=&arrSum(@numOfSpeciesCor);
			    $S2=&arrSum(@numOfReads);

			    if($S0!=0)
			    {
			    	$P0=$X0/$S0;
			    }
			    else
			    {
			    	$P0=0;
			    }

			    if($S1!=0)
			    {
			    	$P1=$X1/$S1;
			    }
			    else
			    {
			    	$P1=0;
			    }

			    if($S2!=0)
			    {
			    	$P2=$X2/$S2;
			    }
			    else
			    {
			    	$P2=0;
			    }





			    splice(@numOfSpecies, 9, 1);
			    splice(@numOfSpeciesCor, 9, 1);
			    splice(@numOfReads, 9, 1);

			    $std0=&standard_deviation(@numOfSpecies);
			    $std1=&standard_deviation(@numOfSpeciesCor);
			    $std2=&standard_deviation(@numOfReads);



			    #to prove that $transPairReadsRef->{9} was deleted successfully
			    my $count0=$#numOfSpecies+1;
			    my $count1=$#numOfSpeciesCor+1;
			    my $count2=scalar (@numOfReads);


			    $m0=&mean(@numOfSpecies);
			    $m1=&mean(@numOfSpeciesCor);
			    $m2=&mean(@numOfReads);

			    if ($std0>0 && $count0>=5) { $Z0=($X0-$m0)/$std0;} else {$Z0=-10;}#by species irrespective of coordinates
			    if ($std1>0 && $count1>=5) { $Z1=($X1-$m1)/$std1;} else {$Z1=-10;}#by species according to coordinates
			    if ($std2>0 && $count2>=5) { $Z2=($X2-$m2)/$std2;} else {$Z2=-10;}

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
		
		
sub SuppComBitSum
{
	my ($transPairSuppSpeciesTotalRef,$transPairSuppReadsTotalRef)=@_;
	my @scaledSuppSpeciesCom=();
	my @scaledSuppReadsCom=();
	my %transPairSuppSpeciesTotal=%{$transPairSuppSpeciesTotalRef};
	my %transPairSuppReadsTotal=%{$transPairSuppReadsTotalRef};
	foreach my $bitValue (keys %transPairSuppSpeciesTotal)
	{
		#print "bits:",$bitValue,"\n";
		#print Dumper $bitValue,"\n";
		my $scaleFactorSpecies=$transPairSuppSpeciesTotal{$bitValue};
		my $scaleFactorReads=$transPairSuppReadsTotal{$bitValue};
		my $i=0;
		#my @bits=split(//,unpack("b*",$bitValue));
		my @bits=split(//,$bitValue);
		#print "scaleFactor:\t",$scaleFactorSpecies,"\t",$scaleFactorReads,"\n";
		foreach my $b ( @bits ) # $bitValue is a list
		{

			#print $i,"\t",$b,"\n";
			$scaledSuppSpeciesCom[$i]+=$b*$scaleFactorSpecies;
			$scaledSuppReadsCom[$i]+=$b*$scaleFactorReads;
			$i++;
			#print $scaledSuppSpeciesCom[$i],"\t",$scaledSuppReadsCom[$i],"\n";
		}

	}

	return \@scaledSuppSpeciesCom,\@scaledSuppReadsCom;
}

sub mean 
{
	my $count=0;
	my(@numbers) =@_;
	foreach (@_) { $count+=$_;}
	return $count/(scalar @_);
}

sub arrSum
{
	my $sum=0;
	my (@numbers) =@_;
	foreach (@_) { $sum+=$_;}
	return $sum;
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
				print "-j  <inputfile2[can be as same as inputfile1]>\n\t";
				print "-o  <outputdir>\n\t";
				print "-b  <bowtie index outputdir>\n\t";
				print "-n  <number of inputfiles>\n\t";
				print "-s  <species name[fly|bombyx|mouse]>\n\t";
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
