#!/bin/sh

#11/05/2013
#WEI WANG
#for Ping-Pong method paper, analyze the data from Yuki lab
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
INDIR=/home/wangw1/data/projects/uava/mouse


OUTDIR=/home/wangw1/data/projects/uava/mouse/masterTablePP
[ ! -f ${OUTDIR} ] && mkdir ${OUTDIR}
LOG=${OUTDIR}/LOG.txt

for i in *v10/*.UA_VA.pp
do
	##
	filename=${i##*/}
	cat $i |grep trans |cut -f1,3,4 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.prefixSpe.txt
	cat $i |grep trans |cut -f1,3,5 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.piSpe.txt
	cat $i |grep trans |cut -f1,3,6 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.pairedReadsSpe.txt
	
	[ ! -s ${OUTDIR}/masterTable.${filename}.prefixSpe.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.prefixSpe.txt ${OUTDIR}/masterTable.${filename}.prefixSpe.txt && rm ${OUTDIR}/${filename}.prefixSpe.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.piSpe.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.piSpe.txt ${OUTDIR}/masterTable.${filename}.piSpe.txt && rm ${OUTDIR}/${filename}.piSpe.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.pairedReadsSpe.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.pairedReadsSpe.txt ${OUTDIR}/masterTable.${filename}.pairedReadsSpe.txt && rm ${OUTDIR}/${filename}.pairedReadsSpe.txt
							
done
