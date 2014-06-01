#!/bin/sh

export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
INDIR=/home/wangw1/data/projects/cd/smRNA/jia_pipeline_results

##cluster annotation (on zlab3)
ANNOA=/home/wangw1/pipeline_dm/common/JB_cluster.42AB.bed
ANNOB=/home/wangw1/pipeline_dm/common/JB_cluster.flam.bed

## bedintersect with cluster ( strand?)

OUTDIR=/home/wangw1/data/projects/cd/smRNA/phasingbyCluster
[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}


parafly_file=${OUTDIR}/intersect.${RANDOM}.para && \
rm -rf $parafly_file

for i in `ls ${INDIR}/*.inserts/*.xkxh.norm.bed.gz`
do 	
	FILE=${i##*/}
	insertsname=`basename $FILE .xkxh.norm.bed.gz`
	inserts=${FILE%%.inserts.*}
	inserts=${inserts}.inserts
	
	candidateflam=${INDIR}/${insertsname}/${insertsname}.xkxh.norm.bed.chrX
	candidate42AB=${INDIR}/${insertsname}/${insertsname}.xkxh.norm.bed.chr2R
	
	mapper42AB=${INDIR}/${insertsname}/${insertsname}.xkxh.norm.bed.42AB
	mapperflam=${INDIR}/${insertsname}/${insertsname}.xkxh.norm.bed.flam
	
	#assume cluster annotation in bed format
	echo "zcat $i |grep -v data |grep chrX |bedtools sort -i stdin > $candidateflam && " >> $parafly_file ; 
	echo "bedtools intersect -a $candidateflam -b $ANNOB -f 0.99 -wa >  ${mapperflam} && rm $candidateflam " >> $parafly_file ; 
	
	echo "zcat $i |grep -v data |grep chr2R |bedtools sort -i stdin > $candidate42AB	&& " >> $parafly_file ;	
	echo "bedtools intersect -a $candidate42AB -b $ANNOA -f 0.99 -wa >  ${mapper42AB} && rm $candidate42AB " >> $parafly_file ;
done
if [ ! -f ${parafly_file}.completed ] || [ -f $parafly_file.failed_commands ]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi
