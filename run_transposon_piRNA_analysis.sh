#!/bin/bash


################
# Major Config #
################
# pipeline version
export smallRNA_downstream_analysis=1.0.0
# this pipeline is still in debug mode
export DEBUG=1 
# pipeline address: if you copy all the files to another directory, this is the place to change; under this directory sits two directories, bin and common. bin stores all the binary executables and common stores all the information of each ORGANISM.
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
# set PATH to be aware of pipeline/bin; this bin directory needs to be searched first
export PATH=${PIPELINE_DIRECTORY}/:$PATH

INDIR=$1 #this is the folder store all pipeline results outmost folders
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log
STEP=1
echo -e "`date "+$ISO_8601"`\tDraw polar histogram of sense fraction" >> $LOG
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.senseFraction ] && \
for i in `ls ${INDIR}/*.inserts/output/*.transposon.list`
do 

	FILE=${i##*/}
	insertsname=`basename $FILE .transposon.list`
	OUTDIR=${INDIR}/transposon_piRNA
	[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}
#polarHistogram for sense fraction
${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/polarHistogram_sense_fraction.r plot_PolarHisto_senseFraction $i ${OUTDIR}
done
[ $? == 0 ] && \	
	touch ${OUT}/.status.${STEP}.transposon_piRNA.senseFraction
STEP=$((STEP+1))

echo -e "`date "+$ISO_8601"`\tDraw length distribution of transposon piRNAs" >> $LOG
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.lendis2 ] && \
OUTDIR=${INDIR}/transposon_piRNA/lendis && \
[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR} && \
paraFile=${OUTDIR}/${RANDOM}.drawlendis2.para && \
for i in `ls ${INDIR}/*.inserts/*.xkxh.transposon.mapper2.gz`
do 	
	FILE=${i##*/}
	insertsname=`basename $FILE .xkxh.transposon.mapper2.gz`
	inserts=${FILE%%.inserts.*}
	inserts=${inserts}.inserts
	nfnnc=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	nfdep=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	
	echo -e "${PIPELINE_DIRECTORY}/lendis2.pl ${i} $OUTDIR nnc $nfnnc &&" >>${paraFile}
	echo -e "RRR ${PIPELINE_DIRECTORY}/R.source plot_lendis2 ${OUTDIR}/$insertsname.xkxh.transposon.mapper2.nnc.lendis2 " >>${paraFile}	
	echo -e "${PIPELINE_DIRECTORY}/lendis2.pl ${i} $OUTDIR seqDep $nfnnc &&" >>${paraFile}
	echo -e "RRR ${PIPELINE_DIRECTORY}/R.source plot_lendis2 ${OUTDIR}/$insertsname.xkxh.transposon.mapper2.seqDep.lendis2 " >>${paraFile}		
done
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA.senseFraction
STEP=$((STEP+1))