#!/bin/bash

export PIPELINE_PATH=/home/wengz/pipelines/smallRNApipeline/pipeline_dm/
export PATH=${PIPELINE_PATH}:${PATH}
ADAPTER=TGGAATTCTCGGGTGCCAAGG
INDIR=$1
OUTDIR=$2
paraFile=${RANDOM}.para
#discard sequences shorter than 18nt
#discard non-adapted sequences
#require minimum adapter length of 8
for fq in `ls *.ovary.fq`
do \
file=${fq##*/}
filename=`basename $file .fq`
#echo -e "fastx_clipper -a $ADAPTER -l 18 -c -M 8 -i $fq -o ${OUTDIR}/${filename}.inserts" >>${paraFile}
echo -ne " Extract_insert_6mer.pl $fq $ADAPTER > ${filename}.insertsout && " >>${paraFile}
echo -ne "inserts2uniqreads.pl ${filename}.insertsout 18 30 > ${filename}.inserts && " >>${paraFile}
echo -e "rm ${filename}.insertsout" >>${paraFile}
done

if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	CPUN=`wc -l $paraFile |cut -f1 -d" "` && \
	ParaFly -c $paraFile -CPU $CPUN -failed_cmds $paraFile.failed_commands
fi
