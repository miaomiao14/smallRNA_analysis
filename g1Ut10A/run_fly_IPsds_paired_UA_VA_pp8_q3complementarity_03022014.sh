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

indexFlag=$1 #this is the folder store all pipeline results outmost folders 
#/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results
fa=/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta

STEP=1

script=${PIPELINE_DIRECTORY}/g1Ut10A/pp8_q2_ww1_zscore_sep_02272014_baseFraction.pl

#run UA_VA analysis
#separate sense and antisense from normbed
#Phil.AubIP.AubCDrescue.unox.ovary.inserts.xkxh.norm.bed.transposons.gz

#uniq bound piRNAs exlude Piwi bound(less piRNA input, to be consistent with the analysis done in 2012)
#use xkxh.norm.bed instead of transposon mappers
INDIR=/home/wangw1/nearline/mpo/IP
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log
declare -a GT=("ago3Hets" "aubHets" "qinHets")
declare -a OX=("unox")
declare -a UNIQ=("uniq")
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans_prefix_complementarity
[ ! -d ${OUTDIR} ] && mkdir $OUTDIR
#touch ${OUT}/.status.${STEP}.all_piRNA_exclude_Piwi.UA_VA
if [ ! -f ${OUT}/.status.${STEP}.all_piRNA_exclude_Piwi.UA_VA_cis_trans_prefix_complementarity ]
then 
	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
			for s in ${UNIQ[@]}
			do
				A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				[ -f ${A} ] && [ -f ${B} ] && \
				for i in {9..25}
				do
					jobname=${t}${s}_${o}_Ago3_Aub.${i}.pp8.q2 && \
					jOUT=${OUTDIR}/${t}${s}_${o}_${i}_Ago3_Aub && \
					[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
					wsize=??
					/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed -a ${fa} -w $wsize -p $i >${jOUT}/${t}${s}_${o}_${i}_Ago3_Aub.pp8.q2.UA_VA_prefix_complementarity.log"
				done
			done
		done
	done

fi

#use xkxh.norm.bed instead of transposon mappers
INDIR=/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans_prefix_complementarity
[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}


echo -e "`date` "+$ISO_8601"\trun pp8 UA and VA..." >> $LOG
#declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "ago3Hets" "aubHets" "qinHets" "nosAgo3CDrescue" "nosAgo3WTrescue")
declare -a GT=("w1")
declare -a OX=("ox" "unox")
declare -a UNIQ=("uniq")

if [ ! -f ${OUT}/.status.${STEP}.all_piRNA.UA_VA_cis_trans_prefix_complementarity ]
then 
	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
			for s in ${UNIQ[@]}
			do
				A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				[ -f ${A} ] && [ -f ${B} ] && \
				for i in {9..25}
				do
				jobname=${t}${s}_${o}_Ago3_Aub.${i}.pp8.q2 && \
				jOUT=${OUTDIR}/${t}${s}_${o}_${i}_Ago3_Aub && \
				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
				wsize=??
				/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed -a ${fa} -w $wsize -p $i >${jOUT}/${t}${s}_${o}_${i}_Ago3_Aub.pp8.q2.UA_VA_prefix_complementarity.log"		
				done
			done
		done
	done
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.all_piRNA.UA_VA_cis_trans_prefix_complementarity
STEP=$((STEP+1))


