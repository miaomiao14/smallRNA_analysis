#!/usr/bin/perl 
use warnings;
use Getopt::Long;
 
=begin comment <runJob.pl notes>
The purpose of this script is to make it very easy to run single use jobs on the 
cluster by generating a one time use script that is submitted after creation. 
 
Options:  
  -i    input command e.g. "muscle -in ensemblSeqs.fa -out ensembl.aln (required)
  -c    number of cores, on hpcc01/02 there are 24 max (required)
  -m    gb's of ram requested, max practical = 45  (required)
  -d    directory to place the temp script in (default = /home/wespisea/tempScripts/
  -t    tag appended to the end of the temp script for identification later (default = "") 
  -test boolean flag, if set script will be crated but qsub command will not be run 
 
 
To adopt this script , please re-name the my $dir to a default
directory with read/write permission and also the email address
 
@2012 Adam Wespiser, University of Massachusetts Medical School
=cut 
 
 
my $test = 0; 
my ($inFile,$dir,$in, $out, $cores, $memory,$tag,$wait,$queue);
GetOptions('i=s' => \$inFile,
           'c=s' => \$cores,
           'm=s' => \$memory,
		   'W=s' => \$wait,
		   'Q=s' => \$queue,
           't=s' => \$tag,
           'd=s' => \$dir,
           'test' => \$test,
           ) or &usage;
 
print $test."\n\n\n"; 
if (!defined $inFile) { print STDERR "\n\-i option was not specified\n";&usage; }
if (!defined $cores)  { print STDERR "\n\-c option was not specified\n";&usage; }
if (!defined $memory) { print STDERR "\n\-m option was not specified\n";&usage; }
if (!defined $wait) { print STDERR "\n\-W option was not specified\n";&usage; }
if (!defined $queue) { print STDERR "\n\-Q option was not specified\n";&usage; }
if (!defined $tag)    { $tag="qsubJob"; } else { $tag = $tag; }
#if (!defined $tag)    { $tag=""; } else { $tag = "_".$tag; }
if (!defined $dir)    { $dir = "/home/aw30w/log/tempScripts"; } else {$dir=~s/\/$//g; }
#if (!defined $dir)    { $dir = "/Users/adam/log/tempScripts"; } else {$dir=~s/\/$//g; }
 
my $abbr = &createAbbr($inFile, $tag);
 
my $heredoc = <<END;
#!/bin/sh
#BSUB -L /bin/bash
#BSUB -n $cores 
#BSUB -R rusage[mem=${memory}] 
#BSUB -q $queue
#BSUB -W $wait
#BSUB -e $dir/${abbr}.err 
#BSUB -o $dir/${abbr}.out
 
$inFile  
END
 
my $file = "$dir/$abbr.temp"; 
 
print "a file has been created :<".$file.">\n"; 
 
open my $fh, '>', $file or die "file cannot be opened\n";
print $fh $heredoc;
close $fh;
print $file."\n"; 
 
if ($test ~~ 0){
  #print "would have run it...\n"; 
  system("bsub < $file"); 
}
 
#use perl's variable hashcode to generate a random hex string 
#(need's to be tested for randomness) 
# another option is to use '$$', perl's var for processID
sub hashifyString{
  my $in = shift;  
  my $f = \$in; 
  $f =~s/(?:SCALAR|\(|\))//g; 
  return $f; 
}
 
sub usage{
  print STDERR "\n\n";
  print STDERR "Usage:\n\n";
  print STDERR "the required argument are:\n";
  print STDERR "perl $0 [-i command to run] [-c number of cores] [-m memory in GB]\n";
  print STDERR "example: perl ~/bin/runJob.pl -c 4 -m 40 -i \"muscle -in .../TOLL/immuneDBFiles/agam_ensHomologs.fasta -out .../TOLL/immuneDBFiles/ensembleResults_afa.al\" -t TOLL_ensembl_take2";
  print STDERR "\n\n\n";
  exit;
}
 
sub createAbbr(){
  my ($abbr,$tag)= @_; 
  $abbr=~s/ //g;
  $abbr=~s/[^a-zA-z]//g;
  if(length($abbr) > 20){
    $abbr = substr $abbr, 0, 15; 
  }
  my @timeData = localtime(time);
  my $timeStamp="_pid".$$;
  #my $timeStamp="_pid".$$."_h=".$timeData[2]."_m=".$timeData[1]; 
  return $tag.$timeStamp; 
}