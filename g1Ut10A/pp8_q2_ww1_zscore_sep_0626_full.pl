#!/usr/bin/perl

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "restrict_digts.pm";
require "Jia.pm";
use File::Basename;
use Compress::Zlib;
#use List::Util qw(sum);
#use PDL;  #it's needed for initialize, it has sum functio inside, clash with the one in Util
#use PDL::Char;
#$PDL::SHARE = $PDL::SHARE; 

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

#06/19/2014
#chrinfo

#06/24/2014
#separate cis basefraction and trans basefraction
#add $NTMT{$l[1]} normalization

#06/26/2014
#suppress ambiguous calculation and output because it runs super slowly
#suppress supplemental pairing check
#report both normalize to target and nonnormalize to target: l[1]

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

 my $supLen=23-$basep; #the length of supplemental region

my %chrSize=();
my $chrFile=$parameters->{chrSize};#"/home/wangw1/pipeline_dm/common/dm3.chromInfo";
open CHR, $chrFile or die "could not find $chrFile: $!";
while(my $line=<CHR>)
{
	chomp $line;
	my @l=split(/\t/,$line);
	$chrSize{$l[0]}=$l[1];
}
close(CHR);

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
			

					
			#store the seq of guide 20nt prefix only; for faster extract the reads number later
			$guidepf{$file}{$dnaseq}+=$reads/$ntm;
			
			#the coordinates of each piRNA species
			$guidepfsplit{$file}{$dnaseq}{$piRNA}{"$chr,$fiveend,$strand"}+=$reads/$ntm; #become 0-based from norm.bed format


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
			
			#store the seq of guide 20nt prefix only; for faster extract the reads number later
			$guidepf{$file}{$dnaseq}+=$reads/$ntm;
			$guidepfsplit{$file}{$dnaseq}{$piRNA}{"$chr,$fiveend,$strand"}+=$reads/$ntm; #become 0-based from norm.bed format

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

		            
		            #store chr, 5'end and strand information separately for each inferred target index
		            my $tstrand="-";
		            $targetpfsplit{$file}{$n}{$str}{$piRNA}{"$chr,$fiveend,$tstrand"}+=$reads/$ntm; #store the strand information for potential guide strand

	            	            
        	}
         	else
         	{
         		if($fileFormat eq "bed")
         		{
		            $start=$bedend-$n-1; #closed
					#$supStart=$start+$basep;
		            $fiveend=$start; #0-based
         		}
         		if($fileFormat eq "normbed")
         		{
         			$start=$bedend-$n-1; #closed
         			#$supStart=$start+$basep;
         			$fiveend=$start;	
         		}
	            
	            #$supStart=$start+$basep;
	            #my $endPostion=$supStart+$supLen;
	            #if ($endPostion <= $chrSize{$chr}) #exclude reads boundary the right boundary
	            #{
		            #my $supStr=substr($genome{$chr},$supStart,$supLen);
		            #if(length($supStr) == $supLen) #exclude boundary reads
		            #{
		            	my $str=substr($genome{$chr},$start,$basep);
		            	$targetpf{$file}{$n}{$str}+=$reads/$ntm;
		            	my $tstrand="+";
		            	#store chr, 5'end and strand information separately for each inferred target index	
		            	$targetpfsplit{$file}{$n}{$str}{$piRNA}{"$chr,$fiveend,$tstrand"}+=$reads/$ntm;
		            #}
	            #}
	            	            
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
	
	
	my %transPairSpecies=();
	my %transPairReads=();
	
	my %cisPairSpeciesNonTNorm=();
	my %cisPairReadsNonTNorm=();
	
	
	my %transPairSpeciesNonTNorm=();
	my %transPairReadsNonTNorm=();
	


	open ZSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.zscore.out";
	open PPSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.pp";
	open PPSEQPAIR, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.ppseq.txt";
	#open PPSEQSUPPVECTOR, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.ppsub.vector.scale.txt";
	
	#open PPUAFRACTION, ">$OUTDIR/$guideStrandFile.$targetStrandFile.$basep.prefix.UA_VA.base.fraction.txt";
	print PPSEQPAIR "pairmode\tguidepiRNAs\tguidepiRNAsReads\ttargetpiRNAs\ttargetpiRNAsReads\n";

	
	foreach ($n=0;$n<$wsize;$n++)
	{
		


		# file1 as ref
		$indexb="$BOUTDIR/$targetStrandFile.$basep.$n";
		$seqFile="$QOUTDIR/$guideStrandFile.$basep.seq";
		$bowtieOut="$MOUTDIR/$guideStrandFile.$targetStrandFile.$basep.$n.bowtie.out";
	 	`[ ! -f $bowtieOut ] && bowtie $indexb -r -a -v 1 -p 8 $seqFile --suppress 1,4,6,7 | grep + > $bowtieOut`;
	   	my %NTMG=();
	   	my %NTMT=();
	   	#my %NTMGP=();
	   	#my %NTMTP=();
	   	open IN, "$MOUTDIR/$guideStrandFile.$targetStrandFile.$basep.$n.bowtie.out";
	   	while(my $line=<IN>)
	   	{
		   	chomp $line;
		   	@l=split(/\t/,$line);
		   	#($strand, $bowtieindex,$seq,$mm)=split(/\t/,$line); 
		    
		    if ($l[3]=~/(\d+):(\w)>(\w)/)
	      	{

		       next if ($1!=0);  # allow 1mm at the 10th position of target strand
		       
		       #$g_0_nt=$3;
		       #$t_9_nt=&revfa($2);
		       
	      	}
	      	#for perfect pair
	      	#else
	      	#{
		    #  	$g_0_nt=substr($l[2],0,1);
		    #  	$t_9_nt=&revfa($g_0_nt);
	      	#}
	      	#$NTMGP{$l[2]}{$g_0_nt.$t_9_nt}++;
	      	$NTMG{$l[2]}++; #normalize guide prefix
		    #$NTMTP{$l[1]}{$g_0_nt.$t_9_nt}++;
		    $NTMT{$l[1]}++; #normalize target inferred index		    		    		    		    		 
	   	}
	   	close(IN);


	   	open IN, "$MOUTDIR/$guideStrandFile.$targetStrandFile.$basep.$n.bowtie.out";
	   	while(my $line=<IN>)
	   	{
	      	chomp $line;
	      	@l=split(/\t/,$line);
	      	my $nnGcorTcor=0;
	      	my $gttotal=0;
	      	
			if ($l[3] eq "" or $l[3]=~/^0:(\w)>(\w)/)
			{
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
		    
				$nnGcorTcor=$nGcor*$nTcor; #total species pairs for this prefix
		      	$gttotal=$nGcorTotalReads*$nTcorTotalReads; #total read pairs for this prefix
			}
	      	if ($l[3] eq "")
	      	{
	      	
	      	   $g_0_nt=substr($l[2],0,1); $t_9_nt=&revfa($g_0_nt);  ##here are different from pp8_q2_ww1.pl
		       #targetpf index; guidepf seq
			   #how many of species start with U?
  		        		       		     		      		       
	       	   	#my %cisFlag=(); #this cisFlag is to mark if there is a cis pair among all the possible pairs(by coordinates) among a paired species	    		       
				foreach my $piQuery (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}} )
			    {
			    	
			    	#my ($piRNAGuide, $piGuideSuppSeq)=split(/,/,$piQuery);

					my $guideQuerySpecies=0;
					$guideQuerySpecies=scalar (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}}); #for this guide piRNA, the number of mapping locus(genomic coordinates) 
					my $guideQueryReads=0;

					map {$guideQueryReads+=$_} values %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}};#for this guide piRNA, the number of reads for all mapping locus


					foreach my $piTargetIndex (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}} ) #indexed piRNA prefix; iterate each piRNA
					{	

						my $targetpiIndexSpecies=0;	
						$targetpiIndexSpecies=scalar (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}});#for this target piRNA, the number of mapping locus(genomic coordinates)
						my $targetpiIndexReads=0;

						map {$targetpiIndexReads+=$_} values %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}}; # the total reads for this target piRNA



						my $cisRecordFlag=0; #flag of cispair
						
						foreach my $Grecord (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}}) #prescan each mapping locus for guide piRNA to find cis target
						{
							if($targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}{$Grecord})
							{
								$cisRecordFlag=1; #if there is any cis pair among multiple genomic loci for a piRNA species pair, then this species pairs is composed of cis pair and ambiguous pair
								last;
							}
						}	
												
						if (! $cisRecordFlag) #if there is no cis pair associated with piQuery, then it can be viewed as a trans pair
			       		{
			       			#trans PingPong pair in species
			       			$transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=1/($NTMG{$l[2]}*$NTMT{$l[1]});
			       			#trans PingPong pair in reads
			       			$transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guideQueryReads*$targetpiIndexReads/($NTMG{$l[2]}*$NTMT{$l[1]});
			       			
			       			#trans PingPong pair in species
			       			$transPairSpeciesNonTNorm{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=1/($NTMG{$l[2]});
			       			#trans PingPong pair in reads
			       			$transPairReadsNonTNorm{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guideQueryReads*$targetpiIndexReads/($NTMG{$l[2]});
			       			
			       			foreach my $Grecord (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}}) #each mapping locus for guide piRNA
							{
					       		

						       		#my $piGuideSuppSeq=$GrecordSup;							
									my $guideQueryCorReads=$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}{$Grecord};
						       		#assuming only one supplemental sequence for each coordinates
						       		my $guideQueryCorReadsNorm=$guideQueryCorReads/$NTMG{$l[2]}; #normalize the guide piRNA reads to number of times mapping from on-the-fly bowtie mapping	
								
									# for each mapping locus for this target piRNA, 
									foreach my $Trecord (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}}) 
									{	
									

											my $targetpiIndexCorReads=$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}{$Trecord} ;#the reads number for this target mapping locus
											my $targetpiIndexCorReadsNorm=	$targetpiIndexCorReads/$NTMT{$l[1]};								
											#no matter is cis or trans, has to do the supplemental sequence comparison
											#my $piTargetIndexSuppSeq=$TrecordSup; #supplemental sequence for target strand
#								       		my $diffstr='';
#											if(length($piTargetIndexSuppSeq) == length($piGuideSuppSeq)) #check if the length of guideSuppSeq is the same as that of the $piTargetIndexSuppSeq
#											{
#							      				my $source = PDL::Char->new($piTargetIndexSuppSeq);
#							      				my $match = PDL::Char->new($piGuideSuppSeq);
#					      				     							
#												my $diff = $match == $source; #bit comparison gives bit results: 0, not eqaul; 1 equal
#												$diffstr=join('',$diff->list);#convert the list to a string
#											}
#											else
#											{
#												for(my $j=0;$j<$supLen;$j++)
#												{
#													$diffstr=$diffstr.'0'; #for default value
#											
#												}
#											}
						       								       		
				      						##print out paired piRNA species      							
				      						print PPSEQPAIR "trans\t$piQuery\t$Grecord\t$guideQueryCorReadsNorm\t$piTargetIndex\t$Trecord\t$targetpiIndexCorReadsNorm\n" if ($n==9);

									}#Trecord

							}#Grecord
			       			
			       		}#if flag=1
			       		else #if there is cis pair for this species pair
			       		{	
							foreach my $Grecord (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}}) #each mapping locus for guide piRNA
							{
					     
					       		
 
						       					
									my $guideQueryCorReads=$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}{$Grecord};

					       		 	#assuming only one supplemental sequence for each coordinates
					       			my $guideQueryCorReadsNorm=$guideQueryCorReads/$NTMG{$l[2]};
					       		
					       			###note: how to define cistargets more accurately?
					       			###in addition to excat match, what if just several nucleotides away?
	
					       			#my $tfiveend=$gfiveend;
					       			#my $tstrand=$gstrand;
								
									# for each mapping locus for this target piRNA, 
									foreach my $Trecord (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}}) 
									{	

											my $targetpiIndexCorReads=$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}{$Trecord} ;#the reads number for this target mapping locus
											my $targetpiIndexCorReadsNorm=	$targetpiIndexCorReads/$NTMT{$l[1]};									
											#no matter is cis or trans, has to do the supplemental sequence comparison
																										       														
											my $corPairReads=$guideQueryCorReads*$targetpiIndexCorReads/($NTMG{$l[2]}*$NTMT{$l[1]});
							       			my $corPairSpecies=1/($NTMG{$l[2]}*$NTMT{$l[1]}*$guideQuerySpecies*$targetpiIndexSpecies);
																		
									       	if($Trecord eq $Grecord) #here it only checks this coordinate associated with $piTargetIndex,if there is a cis pair from target piRNA coordinate pool
									       	{				
	       								       		
									       		$cisPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$corPairSpecies; #cis pair species(in terms of piRNA species) must only have one by coordinate definition; but for different record, it has different cis pair
												$cisPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$corPairReads;	
												
												$cisPairSpeciesNonTNorm{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$corPairSpecies*$NTMT{$l[1]}; #cis pair species(in terms of piRNA species) must only have one by coordinate definition; but for different record, it has different cis pair
												$cisPairReadsNonTNorm{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$corPairReads*$NTMT{$l[1]};	       			      							
				      							##print out paired piRNA species      							
				      							print PPSEQPAIR "cis\t$piQuery\t$Grecord\t$guideQueryCorReadsNorm\t$piTargetIndex\t$Trecord\t$targetpiIndexCorReadsNorm\n" if ($n==9);
				      							print PPSEQPAIR "cisNonTNorm\t$piQuery\t$Grecord\t$guideQueryCorReadsNorm\t$piTargetIndex\t$Trecord\t$targetpiIndexCorReads\n" if ($n==9);
				      							last;		      								      								     							
									       	}#if cis coordinate exists,if Grecord eq Trecord
#									       	else
#								       		{		      							
#				      							#ambiguous pairs if not cis pair
#			      								$ambiguousPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$corPairSpecies;
#			      								$ambiguousPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$corPairReads;
#	      								
#			      								print PPSEQPAIR "ambiguous\t$piQuery\t$Grecord\t$piGuideSuppSeq\t$guideQueryCorReadsNorm\t$piTargetIndex\t$Trecord\t$piTargetIndexSuppSeq\t$targetpiIndexCorReadsNorm\t$diffstr\n" if ($n==9);
#	      								
#			      								$ambiguousPairSuppSpecies{$g_0_nt.$t_9_nt}{$diffstr}+=$corPairSpecies;
#			      								$ambiguousPairSuppReads{$g_0_nt.$t_9_nt}{$diffstr}+=$corPairReads;
#	      								
#			      							} #if Grecord neq Trecord
  	
									} #for each Trecord

							} #for each Grecord

			       		}#else flag=1	

					}#$piTargetIndex
			    }#$piQuery
	      	}#perfect pair

	      	elsif ($l[3]=~/(\d+):(\w)>(\w)/)
	      	{

		       next if ($1!=0);  # allow 1mm at the 10th position of target strand
		       
		       $g_0_nt=$3;
		       $t_9_nt=&revfa($2);
		       	
		       #trans PingPong pair in species
		       $transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$nnGcorTcor/($NTMG{$l[2]}*$NTMT{$l[1]});
		       #trans PingPong pair in reads
		       $transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$gttotal/($NTMG{$l[2]}*$NTMT{$l[1]});
		       	
		       #trans PingPong pair in species
			   $transPairSpeciesNonTNorm{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$nnGcorTcor/($NTMG{$l[2]});
			   #trans PingPong pair in reads
			   $transPairReadsNonTNorm{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$gttotal/($NTMG{$l[2]});
		       	
		       	foreach my $piQuery (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}} )
			    {	 
			    	 
			    	foreach my $piTargetIndex (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}} ) #indexed piRNA prefix; iterate each piRNA
					{    
			       		foreach my $Grecord (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}}) #each mapping locus for guide piRNA
						{
				       		
				
								my $guideQueryCorReads=$guidepfsplit{$guideStrandFile}{$l[2]}{$piQuery}{$Grecord};

				       		 	#assuming only one supplemental sequence for each coordinates
				       			my $guideQueryCorReadsNorm=$guideQueryCorReads/$NTMG{$l[2]};
				       		
							
								# for each mapping locus for this target piRNA, 
								foreach my $Trecord (keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}}) 
								{	

										my $targetpiIndexCorReads=$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{$piTargetIndex}{$Trecord} ;#the reads number for this target mapping locus
										my $targetpiIndexCorReadsNorm=	$targetpiIndexCorReads/$NTMT{$l[1]};									

								
								
			       		
			      						##print out paired piRNA species      							
			      						print PPSEQPAIR "trans\t$piQuery\t$Grecord\t$guideQueryCorReadsNorm\t$piTargetIndex\t$Trecord\t$targetpiIndexCorReadsNorm\n" if ($n==9);
			      						print PPSEQPAIR "transNonTNorm\t$piQuery\t$Grecord\t$guideQueryCorReadsNorm\t$piTargetIndex\t$Trecord\t$targetpiIndexCorReads\n" if ($n==9);


								}#Trecord

						}#Grecord
					}#pitarget
				}#piQuery	
			}#elsif
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
					
#					my $n_of_ambiguousPairSpecies=0;
#					$n_of_ambiguousPairSpecies=scalar (keys %{$ambiguousPairSpecies{$p}{$n}});
#					$n_of_ambiguousPairSpecies= &restrict_num_decimal_digits($n_of_ambiguousPairSpecies,3);
#										    
#					my $n_of_ambiguousPairSpecies_cor=0;
#					map {$n_of_ambiguousPairSpecies_cor+=$_} values %{$ambiguousPairSpecies{$p}{$n}} ;
#					$n_of_ambiguousPairSpecies_cor=&restrict_num_decimal_digits($n_of_ambiguousPairSpecies_cor,3);
#											     
#					my $n_of_ambiguousPairReads=0;
#					map {$n_of_ambiguousPairReads+=$_} values %{$ambiguousPairReads{$p}{$n}} ;
#					$n_of_ambiguousPairReads=&restrict_num_decimal_digits($n_of_ambiguousPairReads,3);
#
#					print PPSCOREUA "$m\tambiguous\t$p\t$n_of_ambiguousPairSpecies\t$n_of_ambiguousPairSpecies_cor\t$n_of_ambiguousPairReads\n";
	
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
			   
			   
			   	#Z-score for individual cispairs, did not normalize to target; 
			   foreach my $p (@matchedpairs)
			   {
			  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%{$cisPairSpeciesNonTNorm{$p}},\%{$cisPairReadsNonTNorm{$p}});
			    #how to normalize $X0{$p}?
			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tcisNonTNorm\t$p\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
			   }

			   #Z-score for individual transpairs; by species, reads and by species irrespective of coordinates,did not normalize to target;	   
			   foreach my $p (@pairs)
			   {
			  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%{$transPairSpeciesNonTNorm{$p}},\%{$transPairReadsNonTNorm{$p}});
			    #how to normalize $X0{$p}?
			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\ttransNonTNorm\t$p\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
			   }
			   
			   
			   #Z-score for individual ambiguouspairs; by species, reads and by species irrespective of coordinates	   
#			   foreach my $p (@matchedpairs)
#			   {
#			  	my ($ZofSpecies,$ZofSpeciesCor,$ZofReads,$PofSpecies,$PofSpeciesCor,$PofReads,$PP10ofSpecies,$PP10ofSpeciesCor,$PP10ofReads,$MofSpecies,$MofSpeciesCor,$MofReads,$StdofSpecies,$StdofSpeciesCor,$StdofReads)=&ZscoreCal(\%{$ambiguousPairSpecies{$p}},\%{$ambiguousPairReads{$p}});
#			    #how to normalize $X0{$p}?
#			    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\tambiguous\t$p\t$wsize\t$basep\t$ZofSpecies\t$ZofSpeciesCor\t$ZofReads\t$PofSpecies\t$PofSpeciesCor\t$PofReads\t$PP10ofSpecies\t$PP10ofSpeciesCor\t$PP10ofReads\t$MofSpecies\t$MofSpeciesCor\t$MofReads\t$StdofSpecies\t$StdofSpeciesCor\t$StdofReads\n"; ##file2 is the guide and file1 is the target
#			   }


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
				



			#close(PPUAFRACTION);
			#close(PPSEQSUPPVECTOR);
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
				print "-c  <chromosome size file>\n\t";
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
						elsif($next_arg eq "-c"){ $parameters->{chrSize}= shift(@ARGV); }

		                else{ print "Invalid argument: $next_arg"; usage(); }
		        }
		}		
