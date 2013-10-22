#!/usr/bin/perl
use File::Basename;
if (@ARGV<1) { print "USAGE: cluster_bucket3_ww.pl <cluster.fa> <n> <norm.bed1> <stats table1> optional: <norm.bed2> <stats table2>\n"; exit 1; }

# build bowtie index from cluster.fa with prefix dmel_cluster
open IN, $ARGV[0]; while(<IN>) { if (/>(.+?)\(/) { $name=$1;}  ##cluster.fa
else { chomp; $size{$name}=length($_);}} ##the length of the cluster, can also be calculated by $end-$start+1

for ($n=1;$n<=$ARGV[1];$n++)
{
  
  $file = fileparse($ARGV[2*$n]);
  $filename{$n}=fileparse($ARGV[2*$n]);
  if ($file=~/(.+).uniqmap.xkxh.norm.bed*/){ $name=$1;} ##sample name
  elsif ($file=~/(.+).xkxh.norm.bed*/){ $name=$1;}  
  #open IN, "/home/xuj1/nearline/datasets/smallRNAs_sequencing/mpo/$name/output/$name"."_stats_table_reads";
  open IN,$ARGV[2*$n+1];
  <IN>;
  while(<IN>) {chomp; split(/\t/); $norm_factor=1000000/$_[3];} ##normlization factor: excluding_ncRNAs (sequencing depth)
  open IN, $ARGV[2*$n]; ##norm.bed
  <IN>;
  $gz = gzopen($ARGV[0], "rb") or die "Cannot open $ARGV[$i]: $gzerrno\n" ;
  while($gz->gzreadline($_) > 0)
  #while (<IN>) 
  { chomp; split(/\t/); $NTM{$_[4]}=$_[6];} ##$_[4] is the reads
  $gz->gzclose();
  `grep -v track $ARGV[2*$n] | cut -f5,6 | uniq.lines+ 0 > $file.uniq.reads`; ##f5,f6 are the reads and their reads count
  `run_bowtie.pl $file.uniq.reads 0 /home/wangw1/pipeline/common/indexes/cluster cluster`;  ##run_bowtie.pl will deal with the bowtie output
  `rm $file.uniq.reads.bowtie.out`;
  `rm $file.uniq.reads`;


%cluster=();
%mapper=();
%seq=();
%len_plus_reads=();
%len_plus_species=();
%len_minus_reads=();
%len_minus_species=();
                   
open IN, "$file.cluster"; ##this is the bowtie output transformed to what we usually work on
while(<IN>)
{
  chomp; split(/\t/);
  $_[2]=~/(.+?)\((.+)\)/;  
  $cluster=$1; $loc{$1}=$2;
  ##$cluster is the number of cluster
  # piRNAs
  if (length($_[0])>=23 && length($_[0])<=29)
  {
    foreach ($_[4]..$_[5])
    {  ##4: start ; 5: end
    $cluster{$cluster}{$_[3]}{$_}+=$_[1]/$NTM{$_[0]}; ###$_[3] + or -  ##$_ is the position
    }
    $mapper{$cluster}.="$_[0]\t$_[1]\tchr0:$_[4]-$_[5]($_[3])\tNA\tNA\tNA\t$NTM{$_[0]}\t1\n";
    $seq{$cluster}{$_[0]}=$_[1];
  }
# only 5'end
#  $cluster{$_[2]}{$_[3]}{$_[4]}+=$_[1]/$NTM{$_[0]} if ($_[3] eq '+');
#  $cluster{$_[2]}{$_[3]}{$_[5]}+=$_[1]/$NTM{$_[0]} if ($_[3] eq '-');

#=========lendis================
  if ($_[3] eq '+' )
  {
  print "$_[0]\n" if (!exists $NTM{$_[0]} || $NTM{$_[0]}==0);
  $len_plus_reads{$cluster}{length($_[0])}+=$_[1]/$NTM{$_[0]};
  $len_plus_species{$cluster}{length($_[0])}+=1/$NTM{$_[0]};
  }
  else
  {
  $len_minus_reads{$cluster}{length($_[0])}+=$_[1]/$NTM{$_[0]};
  $len_minus_species{$cluster}{length($_[0])}+=1/$NTM{$_[0]};
  }
}


#-----------------------------------------------------------------------
#=============seqlogo output==============
 %A=(); %C=(); %G=(); %T=();%count=();
 foreach $cluster (keys %cluster)
 {
    open OUT, "> $file.$cluster.seqlogo";
    if (!exists $seq{$cluster})
    {
       foreach (0..28) { print OUT "0\t";} print OUT "\n";
       foreach (0..28) { print OUT "0\t";} print OUT "\n";
       foreach (0..28) { print OUT "0\t";} print OUT "\n";
       foreach (0..28) { print OUT "0\t";} print OUT "\n";
    }
    else
    {
      foreach $seq (keys %{$seq{$cluster}})
      {
         foreach (0..28)
         {
         $A{$_}+=$seq{$cluster}{$seq} if ( substr($seq,$_,1) eq 'A');
         $C{$_}+=$seq{$cluster}{$seq} if ( substr($seq,$_,1) eq 'C');
         $G{$_}+=$seq{$cluster}{$seq} if ( substr($seq,$_,1) eq 'G');
         $T{$_}+=$seq{$cluster}{$seq} if ( substr($seq,$_,1) eq 'T');
         }
      }
       $i=28;
       foreach (0..28)
        {
          $A{$_}=0 if (!exists $A{$_});$C{$_}=0 if (!exists $C{$_});$G{$_}=0 if (!exists $G{$_});$T{$_}=0 if (!exists $T{$_});
          $count{$_}=$A{$_}+$C{$_}+$G{$_}+$T{$_}; if ($count{$_}==0) {$i=$_-1;last;}
        }
          foreach (0..$i) { print OUT $A{$_}/$count{$_},"\t";} print OUT "\n";
          foreach (0..$i) { print OUT $C{$_}/$count{$_},"\t";} print OUT "\n";
          foreach (0..$i) { print OUT $G{$_}/$count{$_},"\t";} print OUT "\n";
          foreach (0..$i) { print OUT $T{$_}/$count{$_},"\t";} print OUT "\n";
    }
  }

#==========pp output ===============================

foreach $cluster (keys %cluster) {
open OUT,"> $file.$cluster.mapper2";
print OUT $mapper{$cluster};
close(OUT);
$p=0;
$p=`/home/wangw1/git/smallRNA_analysis/pp6_clusterbucket.pl $file.$cluster.mapper2`;
open OUT,"> $file.$cluster.pvalue";
printf OUT "%3f\n",$p;
close(OUT);
}



#======== bucket output===============
@strand=qw /+ -/;
 foreach $c (keys %cluster)
 {
    open OUT1, ">$file.$c";
    open OUT2, ">$file.$c.loc";
    print OUT2 "$loc{$c}\n";
    foreach $strand (@strand)
    {
      foreach $pos (1..$size{$c})
      {
      $cluster{$c}{$strand}{$pos}=0 if (!exists $cluster{$c}{$strand}{$pos});
      print OUT1 "$pos\t",$cluster{$c}{$strand}{$pos}*$norm_factor,"\t$strand\n";
      }
    }
  }


#===========lendis output================================

 foreach $c (keys %cluster) { 
   open OUT1, "> $file.$c.lendis.reads";
   open OUT2, "> $file.$c.lendis.species";
   foreach (18..29) { 
   $len_plus_reads{$c}{$_}=0 if (!exists $len_plus_reads{$c}{$_});
   $len_plus_species{$c}{$_}=0 if (!exists $len_plus_species{$c}{$_});
   $len_minus_reads{$c}{$_}=0 if (!exists $len_minus_reads{$c}{$_});
   $len_minus_species{$c}{$_}=0 if (!exists $len_minus_species{$c}{$_});
   print OUT1 "$_\t",$len_plus_reads{$c}{$_}*$norm_factor,"\t",$len_minus_reads{$c}{$_}*$norm_factor,"\n";
   print OUT2 "$_\t",$len_plus_species{$c}{$_}*$norm_factor,"\t",$len_minus_species{$c}{$_}*$norm_factor,"\n";
   }
}
}

#=============pdf============================================
$N=`grep ">"  $ARGV[0] | wc -l | cut -d ' ' -f1`; 
`RRR /home/wangw1/git/smallRNA_analysis/cluster_bucket_coarse_ww.R plot_cluster_bucket $filename{1} $filename{2} $N`;
#`rm $file.*cluster*`;
