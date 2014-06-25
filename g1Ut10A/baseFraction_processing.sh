#!/bin/sh

export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis

INDIR=$1
OUTDIR=$2

[ ! -f ${OUTDIR} ] && mkdir ${OUTDIR}

for i in ${INDIR}/*.UA_VA.ppseq.txt
do
	filename=${i##*/}
	cat $i |grep trans |cut -f2,3,4,5 |sort -u >${OUTDIR}/${filename%.txt}.guide.
	cat $i |grep trans |cut -f6,7,8,9 |sort -u >
done