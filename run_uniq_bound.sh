#! /usr/bin/env bash
#input: inserts1 inserts2 <inserts3>

export PIPELINEDIR=/home/lees2/pipeline:/home/xuj1/pipeline
A=$1
B=$2
#OUTPUTDIR=$3
echo "Find the shared species between file $A and file $B" >> ${A}.${B}.log
/home/wangw1/bin/inserts_file_methods.py -a $A -b $B -i
echo "calculate the number of shared species between file $A and file $B ..." >>${A}.${B}.log
echo "which is:" >>${A}.${B}.log
echo `wc -l ${A}.sharedA` >>${A}.${B}.log
echo "the number of shared reads in file $A is:" >>${A}.${B}.log
echo `sumcol ${A}.sharedA` >>${A}.${B}.log
echo "the number of shared reads in file $B is:" >>${A}.${B}.log
echo `sumcol ${B}.sharedB` >>${A}.${B}.log
/home/wangw1/bin/inserts_file_methods.py -a $A -b $B -u 
echo "calculate the number of unique species for file $A and file $B ..." >>${A}.${B}.log
echo `wc -l ${A}.uniqA` >>${A}.${B}.log
echo `wc -l ${B}.uniqB` >>${A}.${B}.log
echo "the number of reads in file ${A}.uniqA is:" >>${A}.${B}.log
echo `sumcol ${A}.uniqA` >>${A}.${B}.log
echo "the number of reads in file ${B}.uniqB is:" >>${A}.${B}.log
echo `sumcol ${B}.uniqB` >>${A}.${B}.log

