#!/bin/bash

export PIPELINE_PATH=/home/wengz/pipelines/smallRNApipeline/pipeline_dm/
export PATH=${PIPELINE_PATH}:${PATH}
ADAPTER=TGGAATTCTCGGGTGCCAAGG
INDIR=$1
OUTDIR=$2
paraFile=${RANDOM}.para
######################################
#discard sequences shorter than 18nt##
#discard non-adapted sequences      ##
#require minimum adapter length of 8##
#IT DOES NOT WORK                   ##
#DUE TO THE QUALITY SCORE           ##
######################################
for fq in `ls *.ovary.fq`
do \
file=${fq##*/}
filename=`basename $file .fq`
#echo -e "fastx_clipper -a $ADAPTER -l 18 -c -M 8 -i $fq -o ${OUTDIR}/${filename}.inserts" >>${paraFile}
echo -ne "grep -A 1 \"@\" $fq | grep -v \"@\" | grep -v \"\-\-\"  | uniq.reads+   >${fq}.raw && " >>${paraFile}
echo -ne " Extract_insert_10mer.pl ${fq}.raw $ADAPTER > ${OUTDIR}/${filename}.insertsout && " >>${paraFile}
echo -ne "inserts2uniqreads.pl ${OUTDIR}/${filename}.insertsout 18 30 > ${OUTDIR}/${filename}.inserts && " >>${paraFile}
echo -e "rm ${OUTDIR}/${filename}.insertsout" >>${paraFile}
done

if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	CPUN=`wc -l $paraFile |cut -f1 -d" "` && \
	ParaFly -c $paraFile -CPU $CPUN -failed_cmds $paraFile.failed_commands
fi