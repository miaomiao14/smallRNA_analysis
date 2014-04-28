#!/bin/sh

#11/05/2013
#WEI WANG
#for Ping-Pong method paper, analyze the data from Yuki lab
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
INDIR=/home/wangw1/data/projects/uava/fly


OUTDIR=/home/wangw1/data/projects/uava/fly/masterTablePP
[ ! -f ${OUTDIR} ] && mkdir ${OUTDIR}
LOG=${OUTDIR}/LOG.txt

for i in *v10/*.UA_VA.pp
do
	##
	filename=${i##*/}
	cat $i |grep trans |cut -f1,3,4 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.prefixSpe.txt
	cat $i |grep trans |cut -f1,3,5 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.SpeciesPairs.txt
	cat $i |grep trans |cut -f1,3,6 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.ReadsPairs.txt
		
	cat $i |grep cis |cut -f1,3,4 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.prefixSpe.cis.txt
	cat $i |grep cis |cut -f1,3,5 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.SpeciesPairs.cis.txt
	cat $i |grep cis |cut -f1,3,6 |awk 'BEGIN{OFS="\t"}{print $2,$1,$3}' >${OUTDIR}/${filename}.ReadsPairs.cis.txt
	
	[ ! -s ${OUTDIR}/masterTable.${filename}.prefixSpe.cis.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.prefixSpe.cis.txt ${OUTDIR}/masterTable.drosophila.${filename}.prefixSpePairs.cis.txt && rm ${OUTDIR}/${filename}.prefixSpe.cis.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.SpeciesPairs.cis.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.SpeciesPairs.cis.txt ${OUTDIR}/masterTable.drosophila.${filename}.SpeciesPairs.cis.txt && rm ${OUTDIR}/${filename}.SpeciesPairs.cis.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.ReadsPairs.cis.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.ReadsPairs.cis.txt ${OUTDIR}/masterTable.drosophila.${filename}.ReadsPairs.cis.txt && rm ${OUTDIR}/${filename}.ReadsPairs.cis.txt
	
	
	[ ! -s ${OUTDIR}/masterTable.${filename}.prefixSpe.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.prefixSpe.txt ${OUTDIR}/masterTable.drosophila.${filename}.prefixSpePairs.trans.txt && rm ${OUTDIR}/${filename}.prefixSpe.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.SpeciesPairs.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.SpeciesPairs.txt ${OUTDIR}/masterTable.drosophila.${filename}.SpeciesPairs.trans.txt && rm ${OUTDIR}/${filename}.SpeciesPairs.txt
	[ ! -s ${OUTDIR}/masterTable.${filename}.ReadsPairs.txt ] && ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${filename}.ReadsPairs.txt ${OUTDIR}/masterTable.drosophila.${filename}.ReadsPairs.trans.txt && rm ${OUTDIR}/${filename}.ReadsPairs.txt
							
done
