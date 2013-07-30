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
[ ! -d ${INDIR}/transposon_piRNA ] && mkdir -p ${INDIR}/transposon_piRNA

for i in `ls ${INDIR}/*.inserts/output/*.transposon.list`
do 
	#ln -s $i ${OUTDIR}
	FILE=${i##*/}
	#FILENAME=${FILE%.gz}
	insertsname=`basename $FILE .transposon.list`
	#inserts=${FILE%%.inserts.*}
	#inserts=${inserts}.inserts
	#nfnnc=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	#nfdep=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	#OUTDIR=${INDIR}/transposon_piRNA/${insertsname}
	OUTDIR=${INDIR}/transposon_piRNA
	[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}
	#

#polarHistogram for sense fraction
${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/polarHistogram_sense_fraction.r plot_PolarHisto_senseFraction $i ${OUTDIR}
done