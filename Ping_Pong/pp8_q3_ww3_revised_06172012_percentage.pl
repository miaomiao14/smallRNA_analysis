#!/usr/bin/perl 

BEGIN { unshift @INC,"/home/xuj1/bin";}
require "Statstics.pm";
require "Jia.pm";
use File::Basename;

# what would be the ping-pong z-score be if we require:
# a. "the nucleotides 1-10 of the guide piRNA were perfectly paired with the target piRNA"
# b. "the nucleotides 2-11 of the guide piRNA were perfectly paired with the target piRNA"
# c. "the nucleotides 2-12 of the guide piRNA were perfectly paired with the target piRNA"
# ....
# d. "the nucleotides 2-25 of the guide piRNA were perfectly paired with the target piRNA"
##Jia's norm.bed or mapper2 files are 1-base start and 1-base end!
##In my previous version of pp8 in which I introduce a parameter p
##I should not change the start coordinate for the reference sequences
##In this revised version, I change it, following what Jia did
##06-17-2012

##
open IN, "/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta";
while (<IN>)
{
    if(/>(.+) type/)
    {
        $chr="chr$1";
    }
    else
    {
        chomp;
        $genome{$chr}=$_;
    }
}
$p=$ARGV[3];
#make bowtie index
for($i=0;$i<$ARGV[2];$i++)
{
    my $file=fileparse($ARGV[$i]);
    @namefield=split(/\./,$file);
    $name=$namefield[2]."_".$namefield[1];
    push @argos, $name;
    $file=$name;
    open IN, $ARGV[$i];
    while(<IN>) #complementarity level
    {
        chomp;
        split(/\t/);
	next if (length($_[4])>29 || length($_[4])<23);
	$total{$file}+=$_[5]/$_[6];
        next if (/data/);
            $hash2{$file}{$p}{substr($_[4],0,$p)}+=$_[5]/$_[6];
            for($n=1;$n<=32;$n++) #5'-5'distance
            {
                if($_[3] eq '+')
                {
                    $start=$_[1]+$n-($p+1);
                    $str=substr($genome{$_[0]},$start,$p);
                    $str=&revfa($str);
                    $hash1{$file}{$p}{$n}{$str}+=$_[5]/$_[6];
                }
                else
                {
                    $start=$_[2]-$n; #for bed format file, the end is the actural end coordinate+1 see UCSC bed format
                    $str=substr($genome{$_[0]},$start,$p);
                    $hash1{$file}{$p}{$n}{$str}+=$_[5]/$_[6];
                }
            }
    }


    if($total{$file}>10)
    {
	    open OUT, ">$file.$p.seq";
	    foreach $seq (keys %{$hash2{$file}{$p}})
	    {
		if(length ($seq)==$p)
		{
		    print OUT "$seq\t$hash2{$file}{$p}{$seq}\n";
		}
	    }
	    close(OUT);
	    for ($n=1;$n<=32;$n++)
	    {
		open OUT, ">$file.ref.$p.overlap.$n.distance.fa";
		foreach $ref (keys %{$hash1{$file}{$p}{$n}})
		{
		    if(length ($ref)==$p)
		    {
		    print OUT ">$ref\t$hash1{$file}{$p}{$n}{$ref}\n$ref\n"; #the name of the sequence is the sequence itself
		    }               
		}
		`bowtie-build $file.ref.$p.overlap.$n.distance.fa $file.$p.$n`;
	    }
    }
}

open OUT1,">$p.pp_fraction";
print OUT1 "pp";
print OUT1 map {"\t$_"} @ARGV;
print OUT1 "\n";
open O, ">$p.zscore.out";






# bowtie mapping and score calculating
for ($i=0;$i<$ARGV[2];$i++)
{

    for ($j=0;$j<=$i;$j++)
    {

        $file1=fileparse($ARGV[$i]);
	@namefield=split(/\./,$file1);
	$name1=$namefield[2]."_".$namefield[1];
	$file1=$name1;
        $file2=fileparse($ARGV[$j]);
	@namefield=split(/\./,$file2);
	$name2=$namefield[2]."_".$namefield[1];
	$file2=$name2;

 	#print O "$file1-$file2\t$p";
	#print "$file1-$file2\t$p";
	#print OUT1 "$file1-$file2\t$p";

        if($total{$file1}<10 || $total{$file2}<10)
        {
            print O "$file1-$file2\t$p\t-10\n";
	    print "$file1-$file2\t$p\t-10\n";
            print OUT1 "$file1-$file2\t$p\tNA\n";
        }
        else
        {     
                $X=0;
                $Z=0;
                %score= ();
                $count_N=0;
                open OUT, ">$file1.$file2.$p.pp";
                open OUT2, ">$file1.$file2.$p.ppseq";  #what kind of information does this file store?
                foreach ($n=1;$n<=32;$n++)
                {
                    # file1 as ref, 0 mm
                    `bowtie $file1.$p.$n -r -a -v 1 -p 8 $file2.$p.seq --suppress 1,4,6,7 |grep + > $file1.$file2.$p.$n.bowtie.out`;                    
                    %NTM=undef;
		    open IN, "$file1.$file2.$p.$n.bowtie.out";
                    while (<IN>)
                    {
                        chomp;
                        split(/\t/);
			if ($_[3]=~/(\d+):/) { next if ($1<=9 && $1>=1);}  # no seed mismatches 0 based
			$NTM{$_[2]}++;
		    }
		    close(IN);
		    open IN, "$file1.$file2.$p.$n.bowtie.out";
		    while (<IN>)
                    {
                        chomp;
                        split(/\t/);
			if ($_[3]=~/(\d+):/) { next if ($1<=9 && $1>=1);}  # no seed mismatches 0 based    
		    	
                        $score{$n}+=$hash1{$file1}{$p}{$n}{$_[1]}*$hash2{$file2}{$p}{$_[2]}/$NTM{$_[2]};
                        print OUT2 "NA\t$_[2]\n" if($n==10);  #when 5'-5' distance equals to 10??? $_[2] query sequence
                    }
		    close(IN);
                    # file2 as a ref, 0mm
                    `bowtie $file2.$p.$n -r -a -v 1 -p 8 $file1.$p.seq --suppress 1,4,6,7 |grep + > $file2.$file1.$p.$n.bowtie.out`;
		    %NTM=undef;
                    open IN, "$file2.$file1.$p.$n.bowtie.out";
		    while (<IN>)
                    {
                        chomp;
                        split(/\t/);
			if ($_[3]=~/(\d+):/) { next if ($1<=9 && $1>=1);}  # no seed mismatches 0 based
			$NTM{$_[2]}++;
		    }
		    close(IN);
		    open IN, "$file2.$file1.$p.$n.bowtie.out";
                    while (<IN>)
                    {
                        chomp;
                        split(/\t/);
			if ($_[3]=~/(\d+):/) { next if ($1<=9 && $1>=1);}  # no seed mismatches 0 based
                        $score{$n}+=$hash1{$file2}{$p}{$n}{$_[1]}*$hash2{$file1}{$p}{$_[2]}/$NTM{$_[2]};
                        print OUT2 "$_[2]\tNA\n" if($n==10);                    
                    }
		    close(IN);
                    $score{$n}=0 if (!exists $score{$n});
                    print OUT "$n\t$score{$n}\n";
                    $count_N++ if ($score{$n}>0);
                }
                $X=$score{10};
		%s=();
		for($j=1;$j<=20;$j++)
		{
		   $s{$j}=$score{$j};
		}
		$total_score=0;
		$total_score=&total(values %score);
		if ($total_score!=0)
		{
		$percentage=$X/$total_score;
		}
		else
		{
		   $percentage=0;
		}
		
                delete $score{10};
		$m=&mean(values %score);
                $std=&std(values %score);
                if($std>0 && $count_N>=5)
                {
                    $Z=($X-&mean(values %score))/$std;
                }
                else
                {
                    $Z=-10;
                }
                print O "$file1-$file2\t$p\t$Z\t$percentage\t$m\t$std";
		print "$file1-$file2\t$p\t$Z\t$percentage\t$m\t$std";
		
		$total_score=0;
		$total_score=&total(values %s);
		if ($total_score!=0)
		{
		$percentage=$X/$total_score;
		}
		else
		{
		   $percentage=0;
		}
		delete $s{10};
		$std=&std(values %s);
		$m=&mean(values %s);
		
		if ($X!=0 && $std>0) { $Z=($X-&mean(values %s))/$std;} else {$Z=-10;}
		print "\t$Z\t$percentage\t$m\t$std\t$X\n";
		print O "\t$Z\t$percentage\t$m\t$std\t$X\n";
			     
                if ($Z!=-10)
                {
                    print OUT1 "$file1-$file2\t$p\t", $X/$total{$file1}/$total{$file2}*1000000,"\n"; #what dose this mean?
                }
                else
                {
                    print OUT1 "$file1-$file2\t$p\tNA\n";
                }
                $N1=`/home/xuj1/pipeline/match.pl $file1.$file2.$p.ppseq $file1.$p.seq | sumcol+ 2`;  #$file1.$file2.$p.ppseq OUT2
                chomp($N1);
                `cut -f2 $file1.$file2.$p.ppseq >$file1.$file2.$p.ppseq.temp`;
                $N2=`/home/xuj1/pipeline/match.pl $file1.$file2.$p.ppseq.temp $file2.$p.seq |sumcol+ 2`;
                chomp ($N2);
                open OUT3, ">$file1.$file2.$p.pp_reads_percentage";
                print OUT3 "$N1\t$total{$file1}\t",$N1/$total{$file1},"\t","$N2\t$total{$file2}\t",$N2/$total{$file2},"\n";  
            #}
            
        }
    #print O "\n";
    #print "\n";
    #print OUT1 "\n";       
    }

}

#`rm *.ebwt`;
#`rm *.bowtie.out`;

sub total {
my $count=0;
my(@numbers) =@_;
foreach (@_) { $count+=$_;}
return $count;
}
