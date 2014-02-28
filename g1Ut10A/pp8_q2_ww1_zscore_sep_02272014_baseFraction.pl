#!/usr/bin/perl

BEGIN { unshift @INC,"/home/wangw1/git/smallRNA_analysis/Utils/";}
require "restrict_digts.pm";
require "Jia.pm";
use File::Basename;
use Compress::Zlib;
# rule is p1 doesn't need to pair but p2-20 need
# simplifized version with prefix 20nt
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

#02/27/2014
#fix a bug for background calculation n=17..20, the complementarity actually are 15,14,13,12, can not gurantee 16 nt complementarity
#so for n=1..20, the prefix should be 20 at least  

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




my %totalFirstBase=();
my %totalTenthBase=();
my %totalExpectedTenthBase=();

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
	    if($spe eq "fly")
		{$name=$namefield[2]."_".$namefield[3]."_".$namefield[4]."_".$namefield[11];}
		if($spe eq "bombyx")
	    {$name=$namefield[2]."_".$namefield[12]."_".$namefield[13];}
	    push @argos, $name;
	    $file=$name;
	    
	    &InputFileProcessing($inputfiles[$i],$file);
	     

			$seqFile="$OUTDIR/$file.seq";
			if( ! -s $seqFile )
			{
				open OUT, ">$seqFile";
				foreach my $prefix (keys %{$guidepf{$file}})
				{
					print OUT "$prefix\t$guidepf{$file}{$prefix}\n" if (length($prefix)==20);
				}
			}
			for ($n=0;$n<20;$n++)
			{
				$fa="$OUTDIR/$file.ref.$n.fa";
				$indexb="$OUTDIR/$file.20.$n";
				open OUT, ">$fa";
				foreach my $prefix (keys %{$targetpf{$file}{$n}})
				{
					print OUT ">$prefix\t$targetpf{$file}{$n}{$prefix}\n$prefix\n" if (length($prefix)==20);
				}
				`[ ! -s $indexb ] && bowtie-build $fa $indexb && rm $fa`;
	
			}#for loop of the n

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
#`rm *.fa`;
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


   

			&PingPongProcessing($file2,$file1); #we started from targets and extent the seq to find its reverse complementary guide with certain n of overlap
			if($file1 ne $file2 ) #added on 11/14/2013
			{   				
			   &PingPongProcessing($file1,$file2);    
			}#when file1 and file2 are the same the second iteration is skipped


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
		

		$totalFirstBase{$file}{substr($seq,0,1)}{substr($seq,0,20)}+=$reads/$ntm;
		
		#store the seq of guide 20nt prefix only; for faster extract the reads number later
		$guidepf{$file}{substr($seq,0,20)}+=$reads/$ntm;

		#norm.bed [1,1] -> bed: $2-1,$3
		#bed [0,1) ->norm.bed: $2+1,$3

		#suppose my input is bed; 02/18/2014
		
		#store chr, 5'end and strand information separately for each query 20nt prefix, the populations to find the guide
		my $fiveend=0;	      
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

	  	}
      
      	for (my $n=0;$n<20;$n++)
      	{
      		$totalTenthBase{$file}{$n}{substr($seq,$n,1)}{substr($seq,0,20)}+=$reads/$ntm;
      		my $start=0;
      		my $fiveend=0;
        	if ($strand eq '+') #target strand information
         	{
         		
         		if($fileFormat eq "bed")
         		{
	            	$start=$bedstart+$n+1-20; # $bedstart is 0 based and $start is 0 based, the intermediate end $bedstart+$n+1 is open
	            	$fiveend=$start+20; #open
         		}
         		if($fileFormat eq "normbed")
         		{
         			$start=$bedstart+$n-20;
         			$fiveend=$start+20;#convert to bed, open
         		}
	            my $str=substr($genome{$chr},$start,20); #substr function is 0 based
	            $str=&revfa($str);
	            $targetpf{$file}{$n}{$str}+=$reads/$ntm;
	            

	            
	            $totalExpectedTenthBase{$file}{$n}{substr($str,0,1)}{$str}+=$reads/$ntm;
	            
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
	            my $str=substr($genome{$chr},$start,20);
	            $targetpf{$file}{$n}{$str}+=$reads/$ntm;

	            
	            $totalExpectedTenthBase{$file}{$n}{substr($str,0,1)}{$str}+=$reads/$ntm;
	            
        	}#ifelse
      	}#for
	} #while
	$gz->gzclose();
}

sub PingPongProcessing
{

	my ($guideStrandFile,	$targetStrandFile)=@_;		

	
	open PPUAFRACTION, ">$OUTDIR/$guideStrandFile.$targetStrandFile.UA_VA.base.fraction.txt";
	open PPGSEQ, ">$OUTDIR/$guideStrandFile.ppseq";
	open PPTSEQ, ">$OUTDIR/$targetStrandFile.ppseq";

	
	foreach ($n=0;$n<20;$n++)
	{
		my %pairedFirstBase=();
		my %pairedTenthBase=();
		my %pairedExpectedTenthBase=();
		# file1 as ref
		$indexb="$OUTDIR/$targetStrandFile.20.$n";
		$seqFile="$OUTDIR/$guideStrandFile.seq";
		$bowtieOut="$OUTDIR/$guideStrandFile.$targetStrandFile.$n.bowtie.20.out";
	 	`[ ! -f $bowtieOut ] && bowtie $indexb -r -a -v 1 -p 8 $seqFile --suppress 1,4,6,7 | grep + > $bowtieOut`;
	   	my %NTM=();
	   	open IN, "$OUTDIR/$guideStrandFile.$targetStrandFile.$n.bowtie.20.out";
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
	   	open IN, "$OUTDIR/$guideStrandFile.$targetStrandFile.$n.bowtie.20.out";
	   	while(my $line=<IN>)
	   	{
	      	chomp $line;
	      	@l=split(/\t/,$line);
	      		      	
	      	if ($l[3] eq "")
	      	{
	      	
	      	   	$g_0_nt=substr($l[2],0,1); $t_9_nt=&revfa($g_0_nt);  ##here are different from pp8_q2_ww1.pl

			   
			   	#how many of species start with U?
			   	#$pairedFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]}+=$totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]}/$NTM{$l[2]} ;
		       	$pairedFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]}=$totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]} if($totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]});
			   	$pairedExpectedTenthBase{$targetStrandFile}{$n}{$g_0_nt}{$l[1]}=$totalExpectedTenthBase{$targetStrandFile}{$n}{$g_0_nt}{$l[1]} if($totalExpectedTenthBase{$targetStrandFile}{$n}{$g_0_nt}{$l[1]});#should not be accumulative
			   	# expect to give the same results as %pairedTenthBase
				my $str=&revfa($l[1]);
				$pairedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$str}=$totalTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$str} if($totalTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$str});#should not be accumulative			 

				#$totalTenthBase{$file}{$n}{substr($seq,$n,1)}{substr($seq,0,20)}+=$reads/$ntm;
				print PPGSEQ "$l[2]\n" if ($n==9);
				print PPTSEQ "$str\n" if ($n==9);
	      	}#perfect pair

	      	elsif ($l[3]=~/(\d+):(\w)>(\w)/)
	      	{

		       next if ($1!=0);  # allow 1mm at the 10th position of target strand
		       
		       $g_0_nt=$3;
		       $t_9_nt=&revfa($2);
		      
		       
		       #$pairedFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]}+=$totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]}/$NTM{$l[2]} ;
		       	$pairedFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]}=$totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]} if($totalFirstBase{$guideStrandFile}{$g_0_nt}{$l[2]});
			   	$pairedExpectedTenthBase{$targetStrandFile}{$n}{$g_0_nt}{$l[1]}=$totalExpectedTenthBase{$targetStrandFile}{$n}{$g_0_nt}{$l[1]} if($totalExpectedTenthBase{$targetStrandFile}{$n}{$g_0_nt}{$l[1]});#should not be accumulative
			   	# expect to give the same results as %pairedTenthBase
				my $str=&revfa($l[1]);
				$pairedTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$str}=$totalTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$str} if($totalTenthBase{$targetStrandFile}{$n}{$t_9_nt}{$str});#should not be accumulative			 
				print PPGSEQ "$l[2]\n" if ($n==9);
				print PPTSEQ "$str\n" if ($n==9);
		      
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
	   		
	   		#$totalTenthBase{$file}{substr($seq,0,1)}{substr($seq,0,20)}+=$reads/$ntm;
	   		#total
	   		map {$totalTenthBaseReads{$b}+=$_} values %{$totalTenthBase{$targetStrandFile}{$n}{$b}};
	   		$totalTenthBaseReadsTotal+= $totalTenthBaseReads{$b};
	   		#species
	   		$totalTenthBaseSpecies{$b}=scalar (keys %{$totalTenthBase{$targetStrandFile}{$n}{$b}});
	   		$totalTenthBaseSpeciesTotal+= $totalTenthBaseSpecies{$b};
	   }
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
		
			   #ExpectedTenthBaseFraction
	   my %pairedExpectedTenthBaseReads=();
	   my $pairedExpectedTenthBaseReadsTotal=0;
	   my %pairedExpectedTenthBaseReadsF=();
	   my %pairedExpectedTenthBaseSpecies=();
	   my $pairedExpectedTenthBaseSpeciesTotal=0;
	   my %pairedExpectedTenthBaseSpeciesF=();
	   
	   
	   my %totalExpectedTenthBaseReads=();
	   my $totalExpectedTenthBaseReadsTotal=0;
	   my %totalExpectedTenthBaseReadsF=();
	   my %totalExpectedTenthBaseSpecies=();
	   my $totalExpectedTenthBaseSpeciesTotal=0;
	   my %totalExpectedTenthBaseSpeciesF=();
	   
	   foreach my $b (keys %{$pairedExpectedTenthBase{$targetStrandFile}{$n}})
	   {
	   		map {$pairedExpectedTenthBaseReads{$b}+=$_} values  %{$pairedExpectedTenthBase{$targetStrandFile}{$n}{$b}};
	   		$pairedExpectedTenthBaseReadsTotal+=$pairedExpectedTenthBaseReads{$b};
	   		#species
	   		$pairedExpectedTenthBaseSpecies{$b}=scalar (keys  %{$pairedExpectedTenthBase{$targetStrandFile}{$n}{$b}});
	   		$pairedExpectedTenthBaseSpeciesTotal+=$pairedExpectedTenthBaseSpecies{$b};
	   		
	   		#$totalExpectedTenthBase{$file}{substr($seq,0,1)}{substr($seq,0,20)}+=$reads/$ntm;
	   		#total
	   		map {$totalExpectedTenthBaseReads{$b}+=$_} values %{$totalExpectedTenthBase{$targetStrandFile}{$n}{$b}};
	   		$totalExpectedTenthBaseReadsTotal+= $totalExpectedTenthBaseReads{$b};
	   		#species
	   		$totalExpectedTenthBaseSpecies{$b}=scalar (keys %{$totalExpectedTenthBase{$targetStrandFile}{$n}{$b}});
	   		$totalExpectedTenthBaseSpeciesTotal+= $totalExpectedTenthBaseSpecies{$b};
	   }
	   foreach my $b (keys  %pairedExpectedTenthBaseReads)
	   {			
	   		$pairedExpectedTenthBaseReadsF{$b}=$pairedExpectedTenthBaseReads{$b}/$pairedExpectedTenthBaseReadsTotal;
	   		$pairedExpectedTenthBaseReadsF{$b}=&restrict_num_decimal_digits($pairedExpectedTenthBaseReadsF{$b},4);
	   		$pairedExpectedTenthBaseSpeciesF{$b}=$pairedExpectedTenthBaseSpecies{$b}/$pairedExpectedTenthBaseSpeciesTotal;
	   		$pairedExpectedTenthBaseSpeciesF{$b}=&restrict_num_decimal_digits($pairedExpectedTenthBaseSpeciesF{$b},4);
	   		
	   		$totalExpectedTenthBaseReadsF{$b}=$totalExpectedTenthBaseReads{$b}/$totalExpectedTenthBaseReadsTotal;
	   		$totalExpectedTenthBaseReadsF{$b}=&restrict_num_decimal_digits($totalExpectedTenthBaseReadsF{$b},4);
	   		$totalExpectedTenthBaseSpeciesF{$b}=$totalExpectedTenthBaseSpecies{$b}/$totalExpectedTenthBaseSpeciesTotal;
	   		$totalExpectedTenthBaseSpeciesF{$b}=&restrict_num_decimal_digits($totalExpectedTenthBaseSpeciesF{$b},4);
	   }

		$pairedExpectedTenthBaseSpeciesTotal=&restrict_num_decimal_digits($pairedExpectedTenthBaseSpeciesTotal,4);
		$pairedExpectedTenthBaseReadsTotal=&restrict_num_decimal_digits($pairedExpectedTenthBaseReadsTotal,4);
		

		$totalExpectedTenthBaseSpeciesTotal=&restrict_num_decimal_digits($totalExpectedTenthBaseSpeciesTotal,4);
		$totalExpectedTenthBaseReadsTotal=&restrict_num_decimal_digits($totalExpectedTenthBaseReadsTotal,4);
		
		
		
		
		
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

	    print PPUAFRACTION "$m\ttotal\tt10\t$pairedExpectedTenthBaseSpeciesTotal\t$totalExpectedTenthBaseSpeciesTotal\t$pairedExpectedTenthBaseReadsTotal\t$totalExpectedTenthBaseReadsTotal\n";
		foreach my $b(@bases)
		{
			print PPUAFRACTION "$m\t$b\tt10\t$pairedExpectedTenthBaseSpeciesF{$b}\t$totalExpectedTenthBaseSpeciesF{$b}\t$pairedExpectedTenthBaseReadsF{$b}\t$totalExpectedTenthBaseReadsF{$b}\n";
		}
		
		if($m==10)
		{
			$ppgseq="$OUTDIR/$guideStrandFile.ppseq";
			$seqFile="$OUTDIR/$guideStrandFile.seq";
			$ppgseqm="$OUTDIR/$guideStrandFile.ppseq.reads";
			`match.pl $ppgseq $seqFile >$ppgseqm`;
		
			my @bases=("A","C","G","T");
			#total
			my %totalG1guideStat=();
			%totalG1guideStat=&g1Frac($seqFile);
			my $totalG1Species=0;
			my $totalG1Reads=0;
			my %totalG1SpeciesC=();
			my %totalG1ReadsC=();
			my %totalG1SpeciesF=();
			my %totalG1ReadsF=();
			foreach my $b (keys %totalG1guideStat)
			{
				$totalG1SpeciesC{$b}=scalar (keys %{$totalG1guideStat{$b}});
				$totalG1Species+=$g1Species{$b};
				map {$totalG1ReadsC{$b}+=$_} values %{$totalG1guideStat{$b}};
				$totalG1Reads+=$g1Reads{$b};
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
			%pairedG1guideStat=&g1Frac($ppgseqm);
			
			my $pairedG1Species=0;
			my $pairedG1Reads=0;
			my %pairedG1SpeciesC=();
			my %pairedG1ReadsC=();
			my %pairedG1SpeciesF=();
			my %pairedG1ReadsF=();
			foreach my $b (keys %pairedG1guideStat)
			{
				$pairedG1SpeciesC{$b}=scalar (keys %{$pairedG1guideStat{$b}});
				$pairedG1Species+=$g1Species{$b};
				map {$pairedG1ReadsC{$b}+=$_} values %{$pairedG1guideStat{$b}};
				$pairedG1Reads+=$g1Reads{$b};
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
			
			$pptseq="$OUTDIR/$targetStrandFile.ppseq";
			$seqFile="$OUTDIR/$targetStrandFile.seq";
			$pptseqm="$OUTDIR/$targetStrandFile.ppseq.reads";
			`match.pl $pptseq $seqFile >$pptseqm`;
			
			my %totalT10guideStat=();
			%totalT10guideStat=&t10Frac($seqFile);
			my $totalT10Species=0;
			my $totalT10Reads=0;
			my %totalT10SpeciesC=();
			my %totalT10ReadsC=();
			my %totalT10SpeciesF=();
			my %totalT10ReadsF=();
			foreach my $b (keys %totalT10guideStat)
			{
				$totalT10SpeciesC{$b}=scalar (keys %{$totalT10guideStat{$b}});
				$totalT10Species+=$T10Species{$b};
				map {$totalT10ReadsC{$b}+=$_} values %{$totalT10guideStat{$b}};
				$totalT10Reads+=$T10Reads{$b};
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
			%pairedT10guideStat=&t10Frac($ppgseqm);
			
			my $pairedT10Species=0;
			my $pairedT10Reads=0;
			my %pairedT10SpeciesC=();
			my %pairedT10ReadsC=();
			my %pairedT10SpeciesF=();
			my %pairedT10ReadsF=();
			foreach my $b (keys %pairedT10guideStat)
			{
				$pairedT10SpeciesC{$b}=scalar (keys %{$pairedT10guideStat{$b}});
				$pairedT10Species+=$T10Species{$b};
				map {$pairedT10ReadsC{$b}+=$_} values %{$pairedT10guideStat{$b}};
				$pairedT10Reads+=$T10Reads{$b};
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
		
	     
	     
	}#n=1..20

	
	close(PPTSEQ);
	close(PPGSEQ);
	close(PPUAFRACTION);
	
	
	
	
	
	
	
		
}


sub g1Frac
{
	my $input=$_;
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
	return %g1Stat;
}

sub t10Frac
{
	my $input=$_;
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
	return %t10Stat;
}



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


