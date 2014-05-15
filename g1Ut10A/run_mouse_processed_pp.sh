#!/bin/sh

#11/05/2013
#WEI WANG
#for Ping-Pong method paper, analyze the data from Yuki lab
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
INDIR=/home/wangw1/data/projects/uava/mouse


OUTDIR=/home/wangw1/data/projects/uava/mouse/masterTablePP
[ ! -f ${OUTDIR} ] && mkdir ${OUTDIR}
LOG=${OUTDIR}/LOG.trans.txt

for i in *v10/*.UA_VA.pp
do
	##
	filename=${i##*/}
	cat $i |grep trans |cut -f1,3,4 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.prefixSpe.trans.txt
	cat $i |grep trans |cut -f1,3,5 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.piSpe.trans.txt
	cat $i |grep trans |cut -f1,3,6 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.pairedReadsSpe.trans.txt
	
	[ ! -s ${OUTDIR}/masterTable.${filename}.prefixSpe.trans.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.prefixSpe.trans.txt ${OUTDIR}/masterTable.DC.${filename}.prefixSpePairs.trans.txt && rm ${OUTDIR}/${filename}.prefixSpe.trans.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.piSpe.trans.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.piSpe.trans.txt ${OUTDIR}/masterTable.DC.${filename}.SpeciesPairs.trans.txt && rm ${OUTDIR}/${filename}.piSpe.trans.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.pairedReadsSpe.trans.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.pairedReadsSpe.trans.txt ${OUTDIR}/masterTable.DC.${filename}.ReadsPairs.trans.txt && rm ${OUTDIR}/${filename}.pairedReadsSpe.trans.txt
							
done
