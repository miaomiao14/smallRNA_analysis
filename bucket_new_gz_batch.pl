#!/usr/bin/perl
# Author : Soo Lee (duplexa@gmail.com)

##QUESTIONS## (user input)

#$inputfilename1=&QQQ("What is the name of your 1st input file? Extention must be .transposon.mapper2. (ex. 29SEP08.s1.xkxh.transposon.mapper2)");
#$inputfilename2=&QQQ("What is the name of your 2nd input file? (ex. 29SEP08.s2.xkxh.transposon.mapper2)");
#$inputdir=&QQQ("What is input directory (full path)? (ex. /home/lees2/nearline/29SEP08/)");
#$outputdir=&QQQ("Where should the output files be placed (full path)? (Enter if same as input directory)");
#if($outputdir eq "") { $outputdir = $inputdir; }
#
#$samplename1=&QQQ("Describe your first samplename, in format of italic:plain. (ex. ago3/ago3:  or  ago3:/TM6B   or  :wt)");
#$samplename2=&QQQ("Describe your second samplename likewise.");
#
#$seqdepth1=&QQQ("What is the normalizing factor for your first sample? (for example, sequencing depth or total ncRNA in unit of reads per million) (eg. 1.448)");
#$seqdepth2=&QQQ("What is the normalizing factor for your second sample?");
#
#$email=&QQQ("Your email address");

$inputfilename1=$ARGV[0];
$inputfilename2=$ARGV[1];
$inputdir=$ARGV[2];
$outputdir=$ARGV[3];
if($outputdir eq "") { $outputdir = $inputdir; }
$samplename1=$ARGV[4];
$samplename2=$ARGV[5];

$seqdepth1=$ARGV[6];
$seqdepth2=$ARGV[7];
$email=$ARGV[8];


my ($inputfileprefix1)=$inputfilename1=~/(.+)\.transposon.mapper2.gz$/;
my ($inputfileprefix2)=$inputfilename2=~/(.+)\.transposon.mapper2.gz$/;

$inputdir_exist_message="";
if(! -d $inputdir) {$inputdir_exist_message="(non-existent) *";}
$inputfile1_exist_message="";
if(! -e "$inputdir/$inputfileprefix1.transposon.mapper2.gz") {$inputfile1_exist_message="(non-existent) *";}
$inputfile2_exist_message="";
if(! -e "$inputdir/$inputfileprefix2.transposon.mapper2.gz") {$inputfile2_exist_message="(non-existent) *";}

print "\n\n\nThank you! Check the following information.\n";
print

"Input file 1  $inputfileprefix1.transposon.mapper2.gz $inputfile1_exist_message\n".
"Input file 2  $inputfileprefix2.transposon.mapper2.gz $inputfile2_exist_message\n".
"Input directory   $inputdir $inputdir_exist_message\n".
"Output directory  $outputdir\n".
"Samplename 1    $samplename1\n".
"Samplename 2    $samplename2\n".
"Normalizing factor 1   $seqdepth1\n".
"Normalizing factor 2   $seqdepth2\n".
"Your email address   $email\n\n";

if($inputdir_exist_message ne "" || $inputfile1_exist_message ne "" || $inputfile2_exist_message ne "") {die "Error: Your job cannot be submitted. Please start over.\n\n";}



#$iscorrect=&QQQ("Correct? (Y or N)");
#if($iscorrect ne "Y") {die "Please start over.\n";}




$randomname=""; for(1..30) {$randomname.=(1..9)[rand(10)];}
chomp($home = `echo \$HOME`);
`mkdir -p $home/bucket_sge/`;
$sgefile = "$home/bucket_sge/bucket.$randomname.sge";
print "$sgefile\n";




#transposon list file
$trnlistfilename = "transposon+repeats.list.sorted";
#$trnlistfilename = "transposon.list.lt500+mst40.sorted";
#$trnlistfilename = "transposon+repeats.list.sorted";
#$trnlistfilename = "transposon.list.roo";
#$trnlistfilename = "transposon.list.HMS_Beagle2";
#$trnlistfilename = "transposon.list.1360";

#Make SGE and submit job
#------------------------------------------------

open SGE, ">$sgefile\n";
print SGE
"#!/bin/sh\n".
"\n".
"#\$ -V\n".
"#\$ -pe single 8\n".
"#\$ -l mem_free=31G\n".
"#\$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y\n".
"#\$ -M zzpipeline.admin\@gmail.com\n".
"#\$ -S /bin/bash  \n".
"#\$ -m e \n".
"\n".
"\n".
"#This directory contains the pipeline package\n".
"export PIPELINEDIR1=/home/wengz/pipelines/smallRNApipeline/transposon_bucket\n".
"export PIPELINEDIR2=/home/wengz/pipelines/smallRNApipeline/pipeline_dm\n".
"export PATH=\$PATH:\$PIPELINEDIR1:\$PIPELINEDIR2\n".
"\n".
"#R path\n".
"export PATH=/share/bin/R/bin:/share/bin/R/share:\$PATH\n".
"export LD_LIBRARY_PATH=/share/bin/R/lib64:\$LD_LIBRARY_PATH\n".
"\n".
"#temp, variables\n".
"export BUCKET_INPUTDIR=$inputdir\n".
"export BUCKET_OUTPUTDIR=$outputdir\n".
"export BUCKET_MAPPERPREFIX1=$inputfileprefix1\n".
"export BUCKET_MAPPERPREFIX2=$inputfileprefix2\n".
"export BUCKET_SEQDEPTH1=$seqdepth1\n".
"export BUCKET_SEQDEPTH2=$seqdepth2\n".
"export BUCKET_SAMPLENAME1=$samplename1\n".
"export BUCKET_SAMPLENAME2=$samplename2\n".
"\n".
"export BUCKET_LOGFILE=\$HOME/scratch/jobid_\$JOB_ID/output.log\n".
"export BUCKET_TRNLISTFILE=/home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/$trnlistfilename\n".
"export BUCKET_TRNLENFILE=/home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/transposon_consensus.nodup+mst40+suffix+dodeca.len\n".
"\n".
"\n".
"# Make the directory for the job ID you are running\n".
"# DO NOT DELETE THIS\n".
"mkdir \$HOME/scratch/jobid_\$JOB_ID \n".
"\n".
"#begin\n".
"echo \"Begin\" > \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"\n".
"#COPYING INPUT FROM NEARLINE TO SCRATCH\n".
"echo \"Copying input to scratch\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"cp -r \$BUCKET_INPUTDIR/\$BUCKET_MAPPERPREFIX1.transposon.mapper2.gz \$BUCKET_INPUTDIR/\$BUCKET_MAPPERPREFIX2.transposon.mapper2.gz \$HOME/scratch/jobid_\$JOB_ID/\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"#preparing subdir in outputdir\n".
"\n".
"mkdir \$BUCKET_OUTPUTDIR/bucket_\$JOB_ID\n".
"\n".
"#dividing to FB mappers\n".
"echo \"Dividing mapper to FB mappers\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"FB_gz.pl \$HOME/scratch/jobid_\$JOB_ID/\$BUCKET_MAPPERPREFIX1.transposon.mapper2.gz\n".
"FB_gz.pl \$HOME/scratch/jobid_\$JOB_ID/\$BUCKET_MAPPERPREFIX2.transposon.mapper2.gz\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"#making lendis2\n".
"echo \"Lendis2\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"auto_lendis2.pl  \$BUCKET_MAPPERPREFIX1 \$BUCKET_TRNLISTFILE \$HOME/scratch/jobid_\$JOB_ID/\n".
"auto_lendis2.pl  \$BUCKET_MAPPERPREFIX2 \$BUCKET_TRNLISTFILE \$HOME/scratch/jobid_\$JOB_ID/\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"cp -r \$HOME/scratch/jobid_\$JOB_ID/lendis2 \$BUCKET_OUTPUTDIR/bucket_\$JOB_ID\n".
"\n".
"#Filtering 23-29nt(To make it faster)\n".
"echo \"Filtering 23-29mers\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"auto_filter.pl \$BUCKET_TRNLISTFILE \$BUCKET_MAPPERPREFIX1 \$HOME/scratch/jobid_\$JOB_ID/\n".
"auto_filter.pl \$BUCKET_TRNLISTFILE \$BUCKET_MAPPERPREFIX2 \$HOME/scratch/jobid_\$JOB_ID/\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"#making pairing\n".
"echo \"Pairing(freqsum)\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"auto_pairing.pl \$BUCKET_MAPPERPREFIX1 \$BUCKET_TRNLISTFILE \$HOME/scratch/jobid_\$JOB_ID/ \$BUCKET_SEQDEPTH1\n".
"auto_pairing.pl \$BUCKET_MAPPERPREFIX2 \$BUCKET_TRNLISTFILE \$HOME/scratch/jobid_\$JOB_ID/ \$BUCKET_SEQDEPTH2\n".
"echo \"pp6_bucket\" >> \$BUCKET_LOGFILE\n".
"run \$HOME/scratch/jobid_\$JOB_ID/\$BUCKET_MAPPERPREFIX1.transposon.mapper2 pp6_bucket.pl pp.score\n".
"run \$HOME/scratch/jobid_\$JOB_ID/\$BUCKET_MAPPERPREFIX2.transposon.mapper2 pp6_bucket.pl pp.score\n".
"mv \$HOME/scratch/jobid_\$JOB_ID/*.pp \$HOME/scratch/jobid_\$JOB_ID/*.pp.score \$HOME/scratch/jobid_\$JOB_ID/freqsum\n".
"cd \$HOME/scratch/jobid_\$JOB_ID/freqsum\n".
"cn.sh .pp transposon.mapper2.\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"cp -r \$HOME/scratch/jobid_\$JOB_ID/freqsum \$BUCKET_OUTPUTDIR/bucket_\$JOB_ID\n".
"\n".
"#making seqlogo3\n".
"echo \"Seqlogo3\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"auto_seqlogo3txt \$BUCKET_TRNLISTFILE \$BUCKET_MAPPERPREFIX1 \$HOME/scratch/jobid_\$JOB_ID/ /home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/trneachlenlogo\n".
"auto_seqlogo3txt \$BUCKET_TRNLISTFILE \$BUCKET_MAPPERPREFIX2 \$HOME/scratch/jobid_\$JOB_ID/ /home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/trneachlenlogo\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"#making seqlogo3\n".
"echo \"Seqlogo3 for pingpong pairs\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"auto_seqlogo3txt \$BUCKET_TRNLISTFILE \$BUCKET_MAPPERPREFIX1.pp \$HOME/scratch/jobid_\$JOB_ID/ /home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/trneachlenlogo\n".
"auto_seqlogo3txt \$BUCKET_TRNLISTFILE \$BUCKET_MAPPERPREFIX2.pp \$HOME/scratch/jobid_\$JOB_ID/ /home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/trneachlenlogo\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"cp -r \$HOME/scratch/jobid_\$JOB_ID/fgmatrix \$HOME/scratch/jobid_\$JOB_ID/bgmatrix \$BUCKET_OUTPUTDIR/bucket_\$JOB_ID\n".
"\n".
"\n".
"#consensus mapping\n".
"echo \"Consensus mapping\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"auto_relative_coordination.pl \$BUCKET_TRNLISTFILE \$BUCKET_MAPPERPREFIX1 \$HOME/scratch/jobid_\$JOB_ID/ /home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/trnalnpos\n".
"auto_relative_coordination.pl \$BUCKET_TRNLISTFILE \$BUCKET_MAPPERPREFIX2 \$HOME/scratch/jobid_\$JOB_ID/ /home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/trnalnpos\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"cp -r \$HOME/scratch/jobid_\$JOB_ID/bucket \$BUCKET_OUTPUTDIR/bucket_\$JOB_ID\n".
"\n".
"\n".
"#making seqdepth factor file and samplename file\n".
"echo \"Generating seqdepth and samplename files\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"make_seqdepth_file \$BUCKET_SEQDEPTH1 \$BUCKET_SEQDEPTH2 > \$HOME/scratch/jobid_\$JOB_ID/seqdepthfile\n".
"make_samplename_file \"\$BUCKET_SAMPLENAME1\" \"\$BUCKET_SAMPLENAME2\" > \$HOME/scratch/jobid_\$JOB_ID/samplenamefile\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"#generating pdf\n".
"echo \"Generating pdf\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"\n".
"mkdir -p \$HOME/scratch/jobid_\$JOB_ID/pdf\n".
"RRR /home/wengz/pipelines/smallRNApipeline/pipeline_dm/BUCKETS_new.R totalRNA_bucket \$BUCKET_TRNLISTFILE \$HOME/scratch/jobid_\$JOB_ID/seqdepthfile \$BUCKET_TRNLENFILE \$HOME/scratch/jobid_\$JOB_ID/samplenamefile \$BUCKET_MAPPERPREFIX1 \$BUCKET_MAPPERPREFIX2 \$HOME/scratch/jobid_\$JOB_ID \$BUCKET_MAPPERPREFIX1.\$BUCKET_MAPPERPREFIX2.totalRNA_bucket.pdf\n".
"\n".
"date >> \$BUCKET_LOGFILE\n".
"echo \"done\" >> \$BUCKET_LOGFILE\n".
"echo \"\" >> \$BUCKET_LOGFILE\n".
"\n".
"cp -r \$HOME/scratch/jobid_\$JOB_ID/pdf \$BUCKET_OUTPUTDIR/bucket_\$JOB_ID\n".
"\n".
"\n".
"#END\n".
"echo \"End\" >> \$BUCKET_LOGFILE\n".
"date >> \$BUCKET_LOGFILE\n".
"\n".
"cp \$BUCKET_LOGFILE \$BUCKET_OUTPUTDIR/bucket_\$JOB_ID\n".
"\n".
"mail -s bucket_\$JOB_ID $email < /home/wengz/pipelines/smallRNApipeline/transposon_bucket/commonfiles/bucket_completion_message\n".
"\n";

close SGE;

system("qsub $sgefile");



############################################
## subroutines  ##

sub QQQ {
 my $question = shift @_;
 print "$question\n";
 chomp(my $answer = <STDIN>);
 return $answer;
}


