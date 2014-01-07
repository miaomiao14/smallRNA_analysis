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
BINSIZE=$2 #BINSIZE for cicos plot
OUT=/home/wangw1/isilon_temp/smRNA/transposon_piRNA
LOG=${OUT}/log

echo -e "`date` "+$ISO_8601"\tDraw phasing analysis..." >> $LOG
OUTDIR13=${INDIR}/transposon_piRNA/phasing
[ ! -d $OUTDIR13 ] && mkdir -p ${OUTDIR13}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.phasing ] && \
#paraFile=${OUTDIR13}/${RANDOM}.piRNAphasing.para

for i in `ls ${INDIR}/*.inserts/*.norm.bed.gz`
do
	inputfile=${i##*/}
	samplenamepart=${inputfile#Phil.SRA.*}
	samplename=${samplenamepart%*.xkxh.norm.bed.gz}
	sample=${samplename/ovary.inserts./}
	/home/wangw1/bin/submitsge 8 ${sample} $OUTDIR13 "${PIPELINE_DIRECTORY}/run_distance_analysis.sh -i ${i} -o $OUTDIR13 -t normbed" 
done
touch ${OUT}/.status.${STEP}.transposon_piRNA.phasing