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
OUTDIR1=${INDIR}/transposon_piRNA
echo -e "`date` "+$ISO_8601"\tDraw polar histogram of sense fraction" >> $LOG
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.senseFraction ] && \
[ ! -d $OUTDIR1 ] && mkdir -p ${OUTDIR1} && \
for i in `ls ${INDIR}/*.inserts/output/*.transposon.list`
do 

	FILE=${i##*/}
	insertsname=`basename $FILE .transposon.list`

#polarHistogram for sense fraction
${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/polarHistogram_sense_fraction.r plot_PolarHisto_senseFraction $i ${OUTDIR1}
done
[ $? == 0 ] && \	
	touch ${OUT}/.status.${STEP}.transposon_piRNA.senseFraction
STEP=$((STEP+1))

echo -e "`date` "+$ISO_8601"\tDraw length distribution of transposon piRNAs" >> $LOG
OUTDIR2=${INDIR}/transposon_piRNA/lendis
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.lendis2 ] && \
[ ! -d $OUTDIR2 ] && mkdir -p ${OUTDIR2} && \
paraFile=${OUTDIR2}/${RANDOM}.drawlendis2.para && \
for i in `ls ${INDIR}/*.inserts/*.xkxh.transposon.mapper2.gz`
do 	
	FILE=${i##*/}
	insertsname=`basename $FILE .xkxh.transposon.mapper2.gz`
	inserts=${FILE%%.inserts.*}
	inserts=${inserts}.inserts
	nfnnc=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	nfdep=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	
	echo -ne "${PIPELINE_DIRECTORY}/lendis2.pl ${i} $OUTDIR2 nnc $nfnnc &&" >>${paraFile}
	echo -e "RRR ${PIPELINE_DIRECTORY}/R.source plot_lendis2 ${OUTDIR2}/$insertsname.xkxh.transposon.mapper2.nnc.lendis2 ${insertsname}" >>${paraFile}	
	echo -ne "${PIPELINE_DIRECTORY}/lendis2.pl ${i} $OUTDIR2 seqDep $nfnnc &&" >>${paraFile}
	echo -e "RRR ${PIPELINE_DIRECTORY}/R.source plot_lendis2 ${OUTDIR2}/$insertsname.xkxh.transposon.mapper2.seqDep.lendis2 ${insertsname}" >>${paraFile}		
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA.lendis2
STEP=$((STEP+1))



#declare -a gt_g1_ox=("aubvasAgo3CDrescue.ox" "ago3MutsWW.ox" "aubvasAgo3WTrescue.ox")
#declare -a gt_g1_unox=("aubvasAgo3CDrescue.unox" "ago3MutsWW.unox" "aubvasAgo3WTrescue.unox")
#declare -a gt_g2_ox=("AubCDrescue.ox" "AubMutsWW.ox" "AubWTrescue.ox")
#declare -a gt_g2_unox=("AubCDrescue.unox" "AubMutsWW.unox" "AubWTrescue.unox")
#
#declare -a gt_cor1_ox=("ago3MutsWW.ox" "ago3MutsCJ.ox")
#declare -a gt_cor1_unox=("ago3MutsWW.unox" "ago3MutsCJ.unox")
#
#declare -a gt_cor2_ox=("nosAgo3CDrescue.ox" "aubvasAgo3CDrescue.ox")
#declare -a gt_cor2_unox=("nosAgo3CDrescue.unox" "aubvasAgo3CDrescue.unox")

declare -a GROUPGT=("gt_g1_ox" "gt_g1_unox" "gt_g2_ox" "gt_g2_unox" "gt_cor1_ox" "gt_cor1_unox" "gt_cor2_ox" "gt_cor2_unox")
declare -a gt_g1_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts")
declare -a gt_g1_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts")
declare -a gt_g2_ox=("Phil.SRA.AubCDrescue.ox.ovary.inserts" "Phil.SRA.AubMutsWW.ox.ovary.inserts" "Phil.SRA.AubWTrescue.ox.ovary.inserts")
declare -a gt_g2_unox=("Phil.SRA.AubCDrescue.unox.ovary.inserts" "Phil.SRA.AubMutsWW.unox.ovary.inserts" "Phil.SRA.AubWTrescue.unox.ovary.inserts")

declare -a gt_cor1_ox=("Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.ago3MutsCJ.ox.ovary.inserts")
declare -a gt_cor1_unox=("Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.ago3MutsCJ.unox.ovary.inserts")

declare -a gt_cor2_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts")
declare -a gt_cor2_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts")

echo -e "`date` "+$ISO_8601"\tDraw paired length distribution of transposon piRNAs" >> $LOG
OUTDIR3=${INDIR}/transposon_piRNA/paired_lendis
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.paired.lendis2 ] && \
[ ! -d $OUTDIR3 ] && mkdir -p ${OUTDIR3} && \
paraFile=${OUTDIR3}/${RANDOM}.drawpairedlendis2.para && \
for g in ${GROUPGT[@]}
do
	SUBGROUP=${!g}
	declare -a MAPPER2NNCLENDIS=()
	declare -a MAPPER2UNIQLENDIS=()
	for t in ${SUBGROUP[@]}
	do 
		MAPPER2NNCLENDIS=${MAPPER2NNCLENDIS}" "${OUTDIR2}/${t}.xkxh.transposon.mapper2.nnc.lendis2 
		MAPPER2UNIQLENDIS=${MAPPER2UNIQLENDIS}" "${OUTDIR2}/${t}.uniqmap.xkxh.transposon.mapper2.nnc.lendis2 
	
	done
	#PreMappingList=("${PreMappingList[@]}" "new_target")
	RRR ${PIPELINE_DIRECTORY}/R.source plot_paired_lendis2 ${MAPPER2NNCLENDIS[@]}
	RRR ${PIPELINE_DIRECTORY}/R.source plot_paired_lendis2 ${MAPPER2UNIQLENDIS[@]}
	#echo -e "RRR ${PIPELINE_DIRECTORY}/R.source plot_paired_lendis2 ${MAPPER2NNCLENDIS[@]}" >>${paraFile}	
	#echo -e "RRR ${PIPELINE_DIRECTORY}/R.source plot_paired_lendis2 ${MAPPER2UNIQLENDIS[@]}" >>${paraFile}
done
#[ $? == 0 ] && \
	#ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands &&
	#touch ${OUT}/.status.${STEP}.transposon_piRNA.paired.lendis2
#STEP=$((STEP+1))