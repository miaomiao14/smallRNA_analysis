#!/bin/bash -x


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

script=${PIPELINE_DIRECTORY}/pp6_T_ww_len.pl


INDIR=$1 #this is the folder store all pipeline results outmost folders 
#/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log

STEP=1

declare -a GROUPGT=("AubIP_ago3cdwt_ox" "AubIP_ago3cdwt_unox" \
"AubIP_ago3cdw1_ox" "AubIP_ago3cdw1_unox" \
)

echo -e "`date` "+$ISO_8601"\trun pp8 UA and VA..." >> $LOG
declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "ago3Hets" "aubHets" "qinHets" "nosAgo3CDrescue" "nosAgo3WTrescue")
declare -a OX=("ox" "unox")
declare -a UNIQ=("uniq" "shared")

#uniq
#shared
indexFlag=1 #to indicate we need to build the index or not
OUTDIR=${OUT}/pp6_TOTAL
[ ! -f ${OUTDIR} ] && mkdir -p ${OUTDIR}
#touch ${OUT}/.status.${STEP}.pp6.SRA_vs_SRA
if [ ! -f ${OUT}/.status.${STEP}.pp6.SRA_vs_SRA ] 
then
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do
		for s in ${UNIQ[@]}
		do
			A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.gz
			B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}${s}_${o}_Ago3_Aub.pp6 && \
			jOUT=${OUTDIR}/${t}${s}_${o}_Ago3_Aub && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "${script} ${A} ${B} 2 ${jOUT} >${jOUT}/${t}${s}_${o}_Ago3_Aub.total.pp6.out"
		done
	done
done
#total
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do

			A=${INDIR}/Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.gz
			B=${INDIR}/Phil.SRA.AubIP.${t}.${o}.ovary.inserts/Phil.SRA.AubIP.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}_${o}_Ago3_AubA.pp6 && \
			jOUT=${OUTDIR}/${t}_${o}_Ago3_Aub && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "${script} ${A} ${B} 2 ${jOUT} >${jOUT}/${t}_${o}_Ago3_Aub.total.pp6.out"
	done
done

fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.pp6.SRA_vs_SRA
STEP=$((STEP+1))

