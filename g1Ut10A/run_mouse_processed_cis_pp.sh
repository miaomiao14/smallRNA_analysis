#!/bin/sh

#11/05/2013
#WEI WANG
#for Ping-Pong method paper, analyze the data from Yuki lab
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
INDIR=/home/wangw1/data/projects/uava/mouse


OUTDIR=/home/wangw1/data/projects/uava/mouse/masterTablePP
[ ! -f ${OUTDIR} ] && mkdir ${OUTDIR}
LOG=${OUTDIR}/LOG.cis.txt

for i in *v10/*.UA_VA.pp
do
	##
	filename=${i##*/}
	cat $i |grep cis |cut -f1,3,4 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.prefixSpe.cis.txt
	cat $i |grep cis |cut -f1,3,5 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.piSpe.cis.txt
	cat $i |grep cis |cut -f1,3,6 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.pairedReadsSpe.cis.txt
	
	[ ! -s ${OUTDIR}/masterTable.${filename}.prefixSpe.cis.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.prefixSpe.cis.txt ${OUTDIR}/masterTable.DC.${filename}.prefixSpePairs.cis.txt && rm ${OUTDIR}/${filename}.prefixSpe.cis.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.piSpe.cis.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.piSpe.cis.txt ${OUTDIR}/masterTable.DC.${filename}.SpeciesPairs.cis.txt && rm ${OUTDIR}/${filename}.piSpe.cis.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.pairedReadsSpe.cis.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.pairedReadsSpe.cis.txt ${OUTDIR}/masterTable.DC.${filename}.ReadsPairs.cis.txt && rm ${OUTDIR}/${filename}.pairedReadsSpe.cis.txt
							
done
