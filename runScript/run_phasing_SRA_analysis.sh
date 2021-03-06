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


INDIR=/home/wangw1/data/projects/cd/smRNA/jia_pipeline_results/ #this is the folder store all pipeline results outmost folders
OUT=/home/wangw1/data/projects/cd/smRNA/transposon_piRNA
LOG=${OUT}/log

STEP=1

echo -e "`date` "+$ISO_8601"\tDraw phasing analysis..." >> $LOG
OUTDIR1=${OUT}/phasing
[ ! -d $OUTDIR1 ] && mkdir -p ${OUTDIR1}
if [ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.phasing ] 
then
	for i in `ls ${INDIR}/*.inserts/*.norm.bed.gz`
	do
		inputfile=${i##*/}
		samplenamepart=${inputfile#Phil.SRA.*}
		samplename=${samplenamepart%*.xkxh.norm.bed.gz}
		sample=${samplename/ovary.inserts./}
		/home/wangw1/bin/submitsge 8 ${sample} $OUTDIR1 "${PIPELINE_DIRECTORY}/run_distance_analysis.sh -i ${i} -o $OUTDIR1 -t normbed -r SRA" 
	done
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.transposon_piRNA.phasing
STEP=$((STEP+1))

echo -e "`date` "+$ISO_8601"\tgenerate phasing master table..." >> $LOG
#

OUTDIR2=${OUT}/phasingMaster
[ ! -d $OUTDIR2 ] && mkdir -p ${OUTDIR2}
if [ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.phasing.mastertable ] 
then
	[ -s ${OUTDIR2}/allpiRNAs.allgt.5-5.distance.min.distribution.summary.raw.txt ] && rm ${OUTDIR2}/allpiRNAs.allgt.5-5.distance.min.distribution.summary.raw.txt
	for i in `ls ${OUTDIR1}/*.ovary.inserts.xkxh.norm.bed.gz.5-5.distance.distribution.summary`
	do
		inputfile=${i##*/}
		insertsname=`basename $inputfile .ovary.inserts.xkxh.norm.bed.gz.5-5.distance.distribution.summary`
		samplename=${insertsname#*.SRA.*}
		
	awk -v gt=${samplename} '{OFS="\t"}{print gt,$1,$2}' ${i} >>${OUTDIR2}/allpiRNAs.allgt.5-5.distance.min.distribution.summary.raw.txt
		
		
	done
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR2}/allpiRNAs.allgt.5-5.distance.min.distribution.summary.raw.txt ${OUTDIR2}/allpiRNAs.allgt.5-5.distance.min.distribution.summary.mastertable.txt 
	#transform to fraction (1..100 nt)
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/Phasing/piRNA_distance_plot.r summary_masterTable ${OUTDIR2}/allpiRNAs.allgt.5-5.distance.min.distribution.summary.mastertable.txt ${OUTDIR2}
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.transposon_piRNA.phasing.mastertable
STEP=$((STEP+1))


echo -e "`date` "+$ISO_8601"\tgrouping piRNAs based on their 5'end identity ..." >> $LOG
OUTDIR3=${OUT}/subphasing2
[ ! -d $OUTDIR3 ] && mkdir -p ${OUTDIR3}
if [ ! -f ${OUT}/.status.${STEP}.piRNA.grouping ] 
then
	for i in `ls ${INDIR}/*.inserts/*.norm.bed.gz`
	do
		inputfile=${i##*/}
		samplenamepart=${inputfile#Phil.SRA.*}
		samplename=${samplenamepart%*.xkxh.norm.bed.gz}
		sample=${samplename/ovary.inserts./}
		/home/wangw1/bin/submitsge 2 ${sample} $OUTDIR3 "${PIPELINE_DIRECTORY}/Phasing/piRNAs_grouping.pl ${i} normbed $OUTDIR3 T" 
	done
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.piRNA.grouping 
STEP=$((STEP+1))




echo -e "`date` "+$ISO_8601"\trun phasing of grouped piRNAs based on their 5'end identity ..." >> $LOG
OUTDIR3=${OUT}/subphasing2
[ ! -d $OUTDIR3 ] && mkdir -p ${OUTDIR3}
if [ ! -f ${OUT}/.status.${STEP}.subpiRNA.phasing ] 
then
	for i in `ls ${OUTDIR3}/*.T`
	do
		gzip $i
		inputfile=${i##*/}
		samplenamepart=${inputfile#Phil.SRA.*}
		samplename=${samplenamepart%*.xkxh.norm.bed.*.T}
		sample=${samplename/ovary.inserts./}
	/home/wangw1/bin/submitsge 2 ${sample} $OUTDIR3 "${PIPELINE_DIRECTORY}/run_distance_analysis.sh -i ${i}.gz -o $OUTDIR3 -t normbed -r SRA" 
	done
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.subpiRNA.phasing
STEP=$((STEP+1))

OUTDIR4=${OUT}/phasingMasterT
[ ! -d $OUTDIR4 ] && mkdir -p ${OUTDIR4}
if [ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.T.phasing.mastertable ] 
then
	[ -s ${OUTDIR4}/allpiRNAs.allgt.5-5.distance.min.distribution.summary.raw.txt ] && rm ${OUTDIR4}/allpiRNAs.allgt.5-5.distance.min.distribution.summary.raw.txt
	for i in `ls ${OUTDIR1}/*.ovary.inserts.xkxh.norm.bed.T.gz.5-5.distance.distribution.summary`
	do
		inputfile=${i##*/}
		insertsname=`basename $inputfile .ovary.inserts.xkxh.norm.bed.T.gz.5-5.distance.distribution.summary`
		samplename=${insertsname#*.SRA.*}
		
	awk -v gt=${samplename}.T '{OFS="\t"}{print gt,$1,$2}' ${i} >>${OUTDIR4}/allpiRNAs.allgt.T.5-5.distance.min.distribution.summary.raw.txt
		
		
	done
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR4}/allpiRNAs.allgt.T.5-5.distance.min.distribution.summary.raw.txt ${OUTDIR4}/allpiRNAs.allgt.T.5-5.distance.min.distribution.summary.mastertable.txt 
	#transform to fraction (1..100 nt)
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/Phasing/piRNA_distance_plot.r summary_masterTable ${OUTDIR4}/allpiRNAs.allgt.T.5-5.distance.min.distribution.summary.mastertable.txt ${OUTDIR4}
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.transposon_piRNA.T.phasing.mastertable
STEP=$((STEP+1))




#/home/hanb/scratch/cd/smallRNA_pipeline_output/Phil.SRA.ago3MutsAubMuts.ox.ovary/intersect_piRNA_length/*.insert
#/home/hanb/scratch/cd/smallRNA_pipeline_output/Phil.SRA.w1.ox.ovary/intersect_piRNA_length/*.insert
#seqlogo_ww $i 