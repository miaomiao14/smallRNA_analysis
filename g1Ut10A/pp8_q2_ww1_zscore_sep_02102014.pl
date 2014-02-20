#!/usr/bin/perl
BEGIN { unshift @INC,"/home/xuj1/bin/";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;
use Compress::Zlib;
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
my $indexFlag=$parameters->{indexflag};
my $fileFormat=$parameters->{format};

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



my %targetpf=();
my %guidepf=(); 

my %targetpfsplit=();
my %guidepfsplit=();

my %total=();
#main, preprocessing
if($indexFlag)
{
	# make bowtie index
	for ($i=0; $i<$numOfInput; $i++)
	{
		my $file=fileparse($inputfiles[$i]);
	    @namefield=split(/\./,$file);
	    if($spe eq "fly")
		{$name=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
	    push @argos, $name;
	    $file=$name;
	    
	    &InputFileProcessing($inputfiles[$i],$file);
	     
		if ($total{$file}>10)
		{
			$seqFile="$OUTDIR/$file.seq";
			if( ! -s $seqFile )
			{
				open OUT, ">$seqFile";
				foreach my $prefix (keys %{$guidepf{$file}})
				{
					print OUT "$prefix\t$guidepf{$file}{$prefix}\n" if (length($prefix)==16);
				}
			}
			for ($n=0;$n<20;$n++)
			{
				$fa="$OUTDIR/$file.ref.$n.fa";
				$indexb="$OUTDIR/$file.$n";
				open OUT, ">$fa";
				foreach my $prefix (keys %{$targetpf{$file}{$n}})
				{
					print OUT ">$prefix\t$targetpf{$file}{$n}{$prefix}\n$prefix\n" if (length($prefix)==16);
				}
				`bowtie-build $fa $indexb && rm $fa`;
	
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
    	if($spe eq "fly")
		{$name=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
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
   	    if($spe eq "fly")
		{$name1=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name1=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
		$file1=$name1;
		$file2=fileparse($inputfiles[$j]);
		@namefield=split(/\./,$file2);
   	    if($spe eq "fly")
		{$name2=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name2=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
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
		
		next if (length($seq)>29 || length($seq)<23);
		next if (/data/);
		
		$total{$file}+=$reads/$ntm; #no nta in norm.bed format
		
		#store the seq of guide 16nt prefix only; for faster extract the reads number later
		$guidepf{$file}{substr($seq,0,16)}+=$reads/$ntm;

		#norm.bed [1,1] -> bed: $2-1,$3
		#bed [0,1) ->norm.bed: $2+1,$3

		#suppose my input is bed; 02/18/2014
		
		#store chr, 5'end and strand information separately for each query 16nt prefix, the populations to find the guide
		my $fiveend=0;	      
		if($strand eq '+')
		{	
			
			if($fileFormat eq "bed")
			{ 
				$fiveend=$bedstart; #0-based,closed
			}
			if($fileFormat eq "normbed")
			{
				$fiveend=$bedstart;#1-based,closed
			}
			
			$guidepfsplit{$file}{substr($seq,0,16)}{"$chr,$fiveend,$strand"}+=$reads/$ntm; #become 0-based from norm.bed format
			
			#if store the start of all query
			#my $queryStart=$bedstart-1;
	  	}
	  	else
	  	{
	  		if($fileFormat eq "bed")
			{
	  			$fiveend=$bedend; #open
			}
			if($fileFormat eq "normbed")
			{
				$fiveend=$bedend;#1-based,closed
			}
	  		$guidepfsplit{$file}{substr($seq,0,16)}{"$chr,$fiveend,$strand"}+=$reads/$ntm; #become 0-based from norm.bed format
	  		#my $queryStart=$bedstart-1;
	  	}
      
      	for (my $n=0;$n<20;$n++)
      	{
      		my $start=0;
      		my $fiveend=0;
        	if ($strand eq '+') #target strand information
         	{
         		
         		if($fileFormat eq "bed")
         		{
	            	$start=$bedstart+$n+1-16; # $bedstart is 0 based and $start is 0 based, the intermediate end $bedstart+$n+1 is open
	            	$fiveend=$start+16; #open
         		}
         		if($fileFormat eq "normbed")
         		{
         			$start=$bedstart+$n-16;
         			$fiveend=$start+16;#closed?
         		}
	            my $str=substr($genome{$chr},$start,16); #substr function is 0 based
	            $str=&revfa($str);
	            $targetpf{$file}{$n}{$str}+=$reads/$ntm;
	            
	            #store chr, 5'end and strand information separately for each guide 16nt prefix
	            my $tstrand="-";
	            $targetpfsplit{$file}{$n}{$str}{"$chr,$fiveend,$tstrand"}+=$reads/$ntm; #store the strand information for guide strand
	            #my $indexStart=$start;
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
	            my $str=substr($genome{$chr},$start,16);
	            $targetpf{$file}{$n}{$str}+=$reads/$ntm;
	            my $tstrand="+";
	            #store chr, 5'end and strand information separately for each guide 16nt prefix	
	            $targetpfsplit{$file}{$n}{$str}{"$chr,$fiveend,$tstrand"}+=$reads/$ntm;
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
	my %X0=(); 
	my %Z0=(); 
	my %allPairReads=(); 
	my %suv=(); 
	my %count_N0=();
	my %species=(); 
	my %speciesn10=();
	my %cisPairSpecies=(); my %cisPairReads=(); my %transPairSpecies=(); my %transPairReads=();my %cisPair10Species=();my %transPair10Species=();
	
	
	open ZSCORE, ">$OUTDIR/$guideStrandFile.$targetStrandFile.zscore.out";
	open ZSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.UA_VA.zscore.out";
	 
	open PPSCORE, ">$OUTDIR/$guideStrandFile.$targetStrandFile.pp";
	open PPSCOREUA, ">$OUTDIR/$guideStrandFile.$targetStrandFile.UA_VA.pp";
	

	open PPSEQ, ">$OUTDIR/$guideStrandFile.$targetStrandFile.ppseq";
	
	foreach ($n=0;$n<20;$n++)
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
	      	
	      	$guidetotal=$guidepf{$guideStrandFile}{$l[2]};
		    $targettotal=$targetpf{$targetStrandFile}{$n}{$l[1]};
	      
	      	if ($l[3] eq "")
	      	{
	      	
	      	   $g_0_nt=substr($l[2],0,1); $t_9_nt=&revfa($g_0_nt);  ##here are different from pp8_q2_ww1.pl
		       #targetpf index; guidepf seq
		       
		       
		       
		       $allPairReads{$g_0_nt.$t_9_nt}{$n}+=$targettotal*$guidetotal/$NTM{$l[2]};
		       $score{$n}+=$targettotal*$guidetotal/$NTM{$l[2]};
		       
		       $speciesn10{$g_0_nt.$t_9_nt}{$l[2]}=1 if ($n==9); ###
		       $species{$g_0_nt.$t_9_nt}{$n}{$l[2]}=1 ; #this was wrong, has to add {$n}, otherwise accumulative
		       #the sum of $cisPairSpecies and $transPairSpecies irrespective of coordinates
		       
		       
	      		
		       #separate trans-targets from cis-targets for perfect g1t10 pair  	
		       
		       #$guidepfsplit{$file}{substr($seq,0,16)}{$chr,$fiveend,$strand}+=$reads/$ntm; #become 0-based from 1-based norm.bed format
		       #$targetpfsplit{$file}{$n}{$str}{$chr,$fiveend,"+"}+=$reads/$ntm;#$fiveend is 0-based
		       
		       #how to sort guidepfsplit and targetpfsplit
		       #no need to sort if use hash key
		       
		       
		       foreach my $record (keys %{$guidepfsplit{$guideStrandFile}{$l[2]}} )
		       {
		       		my ($chr,$gfiveend,$gstrand)=split(/,/,$record);
		       		
		       		
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
		       		
		       		
		       		
		       		if($targetpfsplit{$targetStrandFile}{$n}{$l[1]}{"$chr,$tfiveend,$tstrand"})
		       		{
		       			$cisPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targetpfsplit{$targetStrandFile}{$n}{$l[1]}{"$chr,$tfiveend,$tstrand"}/$NTM{$l[2]};
		       			
		       			$cisPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}=1 ; #cis pair species must only have one by coordinate definition
		       			
		       			#with the same guide 16 nt prefix,there might be multiple trans-targets with 16nt complementarity, (originally it was viewed only one)
		       			#trans PingPong pair in species
		       			$transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}})-1 ;
		       			$transPair10Species{$g_0_nt.$t_9_nt}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}})-1 if($n==9) ;
		       			#trans PingPong pair in reads
		       			$transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*($targettotal - $targetpfsplit{$targetStrandFile}{$n}{$l[1]}{"$chr,$tfiveend,$tstrand"})/$NTM{$l[2]};
		       			
		       		}
		       		else
		       		{
		       			#trans PingPong pair in species
		       			$transPairSpecies{$g_0_nt.$t_9_nt}{$n}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}});
		       			#trans PingPong pair in reads
		       			$transPair10Species{$g_0_nt.$t_9_nt}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}}) if($n==9) ;
		       			$transPairReads{$g_0_nt.$t_9_nt}{$n}{$l[2]}+=$guidepfsplit{$guideStrandFile}{$l[2]}{$record}*$targettotal/$NTM{$l[2]};
		       		}

		       }
		       					       
		       #store target seq from query populations, as it's from bowtie output, by default, it has a partner
		       print PPSEQ "$l[2]\n" if ($n==9);
	      	}#perfect pair

	      	elsif ($l[3]=~/(\d+):(\w)>(\w)/)
	      	{

		       next if ($1!=0);  # allow 1mm at the 10th position of target strand
		       $t_9_nt=&revfa($2);
		       
		       $allPairReads{$3.$t_9_nt}{$n}+=$targettotal*$guidetotal/$NTM{$l[2]};
		       $score{$n}+=$targettotal*$guidetotal/$NTM{$l[2]}; ###species of pairs taking the varied coordinates
		       
		       $speciesn10{$3.$t_9_nt}{$l[2]}=1 if ($n==9); ###species of seq pairs, not count different coordinates
		       $species{$3.$t_9_nt}{$n}{$l[2]}=1 ;
		       
		       #trans PingPong pair in species
		       $transPairSpecies{$3.$t_9_nt}{$n}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}});
		       $transPair10Species{$3.$t_9_nt}{$l[2]}=scalar(keys %{$targetpfsplit{$targetStrandFile}{$n}{$l[1]}}) if($n==9);
		       #trans PingPong pair in reads
		       $transPairReads{$3.$t_9_nt}{$n}{$l[2]}+=$guidetotal*$targettotal/$NTM{$l[2]};
		       
		       print PPSEQ "$l[2]\n" if ($n==9);
			}
		} #while
		close(IN);
	   
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
			$n_of_species=scalar (keys %{$species{$p}{$n}});
		     
			my $n_of_cisPairSpecies=0;
			$n_of_cisPairSpecies=scalar (keys %{$cisPairSpecies{$p}{$n}});					    
			my $n_of_cisPairSpecies_cor=0;
			map {$n_of_cisPairSpecies_cor+=$_} values %{$cisPairSpecies{$p}{$n}} ;					     
			my $n_of_cisPairReads=0;
			map {$n_of_cisPairReads+=$_} values %{$cisPairReads{$p}{$n}} ;

			print PPSCOREUA "$m\tcis\t$p\t$n_of_cisPairSpecies\t$n_of_cisPairSpecies_cor\t$n_of_cisPairReads\t$n_of_species\t$allPairReads{$p}{$n}\n";

	   }
	     #for all pairs, trans only 
		foreach my $p (@pairs)
		{
		     $allPairReads{$p}{$n}=0 if (!exists $allPairReads{$p}{$n});
		     $n_of_species=scalar (keys %{$species{$p}{$n}});
		     
		     
		     my $n_of_transPairSpecies=0;					     					   
		     $n_of_transPairSpecies=scalar (keys %{$transPairSpecies{$p}{$n}});
		     my $n_of_transPairSpecies_cor=0;
		     map {$n_of_transPairSpecies_cor+=$_} values %{$transPairSpecies{$p}{$n}};					     
		     my $n_of_transPairReads=0;
		     map {$n_of_transPairReads+=$_} values %{$transPairReads{$p}{$n}} ;
		     
		     print PPSCOREUA "$m\ttrans\t$p\t$n_of_transPairSpecies\t$n_of_transPairSpecies_cor\t$n_of_transPairReads\t$n_of_species\t$allPairReads{$p}{$n}\n";
		     
		     $count_N0{$p}++ if ($n_of_transPairSpecies>0);
	     }
	     
	     
	     
	     
	}#n=1..20
	   $X=$score{9}; delete $score{9};
	   $std=&std(values %score); 
	   if ($std>0 && $count_N>=5) { $Z=($X-&mean(values %score))/$std;} else {$Z=-10;}
	   print ZSCORE "$guideStrandFile\-$targetStrandFile\t$Z\t";
	   print "$guideStrandFile\-$targetStrandFile\t$Z\t";
	   
	   #Z-score for all transpairs; by species irrespective of coordinates
	   
	   foreach my $p (@pairs)
	   {
	    $X0{$p}=scalar (keys %{$transPairSpecies{$p}{9}}); #by species irrespective of coordinates
	    map {$X1{$p}+=$_} (values %{$transPairSpecies{$p}{9}}); #by species according to coordinates
	    map {$X2{$p}+=$_} (values %{$transPairReads{$p}{9}});  #Z-score for all transpairs; by reads
	    
	    
	    $n_of_species=scalar (keys %{$transPair10Species{$p}}); #total number of species irrespective of coordinates
	    my $n_of_species_cor=0;
	    map {$n_of_species_cor+=$_} (values %{$transPair10Species{$p}});#total number of species respective of coordinates
	    my $n_of_reads=0;
	    map {$n_of_reads+=$_} (values %{$transPairReads{$p}{9}}); ###why it is not equal to $X2{$p}? pls calculate it before it's deleted!!
	    
	    delete $transPairSpecies{$p}{9}; ##???
	    delete $transPairReads{$p}{9};###???
	    
	    my @numOfSpecies=();
	    my @numOfSpeciesCor=();
	    my @numOfReads=();
	    
	    for(my $i=1;$i<=20;$i++)
	    {
	    	my $n1=scalar (keys %{$transPairSpecies{$p}{$i}});
	    	my $n2=scalar (keys %{$transPairReads{$p}{$i}});
	    	
	    	if($n1>0) #$transPairSpecies{$p}{9} was deleted
	    	{
	    		push @numOfSpecies, scalar (keys %{$transPairSpecies{$p}{$i}}) ;
	    		my $XnSpecies=0;
	    		map { $XnSpecies+=$_ } (values %{$transPairSpecies{$p}{$i}});
	    		push @numOfSpeciesCor, $XnSpecies ;
	    	}
			if($n2>0) #$transPairReads{$p}{9} was deleted
			{
				my $XnReads=0;
				map { $XnReads+=$_ } (values %{$transPairReads{$p}{$i}}); 
				push @numOfReads, $XnReads ;
			}				    	
	    }
	    
	    $std0=&std(@numOfSpecies);
	    $std1=&std(@numOfSpeciesCor);
	    $std2=&std(@numOfReads);
	    
	    #to prove that $transPairReads{$p}{9} was deleted successfully
	    #my $temp1=$#numOfSpecies+1;
	    #my $temp2=$#numOfSpeciesCor+1;
	    #my $temp3=$#numOfReads+1;
	    
	    $m0=&mean(@numOfSpecies);
	    $m1=&mean(@numOfSpeciesCor);
	    $m2=&mean(@numOfReads);
	    
	    if ($std0>0 && $count_N0{$p}>=5) { $Z0{$p}=($X0{$p}-$m0)/$std0;} else {$Z0{$p}=-10;}#by species irrespective of coordinates
	    if ($std1>0 && $count_N0{$p}>=5) { $Z1{$p}=($X1{$p}-$m1)/$std1;} else {$Z1{$p}=-10;}#by species according to coordinates
	    if ($std2>0 && $count_N0{$p}>=5) { $Z2{$p}=($X2{$p}-$m2)/$std2;} else {$Z2{$p}=-10;}
	    
	    
	   
	    
	    #how to normalize $X0{$p}?
	    print ZSCOREUA "$guideStrandFile\-$targetStrandFile\t$p\t$Z0{$p}\t$Z1{$p}\t$Z2{$p}\t$X0{$p}\t$X1{$p}\t$X2{$p}\t$n_of_species\t$n_of_species_cor\t$n_of_reads\t$m0\t$std0\t$m1\t$std1\t$m2\t$std2\n"; ##file2 is the guide and file1 is the target
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

sub usage
{
        print "\nUsage:$0\n\n\t";
        print "REQUIRED\n\t";
        print "inputfile1 inputfile2 n type(fly|bombyx) outdir indexflag\n\t";
        print "-i  <inputfile1>\n\t";
		print "-j  <inputfile2[can be as same as inputfile1]>\n\t";
		print "-o  <outputdir>\n\t";
		print "-n  <number of inputfiles>\n\t";
		print "-s  <species name[fly|bombyx]>\n\t";
		print "-d  <flag of index[0|1]>\n\t";
		print "-f  <flag of index[bed|normbed]>\n\t";
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
                elsif($next_arg eq "-d"){ $parameters->{indexflag} = shift(@ARGV); }
                elsif($next_arg eq "-f"){ $parameters->{format} = shift(@ARGV); }

                else{ print "Invalid argument: $next_arg"; usage(); }
        }
}

