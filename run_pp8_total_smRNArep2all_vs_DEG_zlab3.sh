#!/bin/bash -x
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp8_ww_smRNA_vs_DEG.pl
INDIR=/home/wangw1/data/projects/cd/smRNA/pp8_SRADEG_datasets
smRNAINDIR=/home/wangw1/data/projects/cd/smRNA/pp8_SRADEG_datasets_rep2
degraINDIR=${INDIR}/degradome/pipeline_output_12262013
OUT0=/home/wangw1/data/projects/cd/smRNA/pp8_smRNAallrep2_vs_degradometrnoutcluster_total_06052014

#note: use window size 16 for complementarity 16 (updated into script)

#use all small RNAs 
declare -a FEATURE=("FLY_TRN_ALL" "FLY_TRN_ALL_IN_CLUSTER" "FLY_TRN_ALL_OUT_CLUSTER")
#g=$5
#c=$6 #cpu
declare -a GROUPGT=("ago3Mutsrep2" "aubvasAgo3CDrescuerep2" "aubvasAgo3WTrescuerep2" "aubMutsrep2" "AubCDrescuerep2" "AubWTrescuerep2")
#declare -a GROUPGT=("ago3MutsWW" "ago3MutsAubMuts" )
#no WW, no rep2

STEP=1

[ ! -d ${OUT0} ] && mkdir -p ${OUT0}
#step 1
#OUT=${OUT0}
#touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8
if [ ! -f ${OUT0}/.status.${STEP}.SRA_DEG.pp8 ] 
then
	for g in "${GROUPGT[@]}"
	do
	if [ $g == "ago3Mutsrep2" ] || [ $g == "aubMutsrep2" ]
	then 
		gt=${g%rep2}WW
	else 
		gt=${g%rep2}
	
	fi
	for f in "${FEATURE[@]}"
	do
		OUTDIR=${OUT0}/${g}_${f}
		[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
		
		LOG=${OUTDIR}/${g}.log	
		#Phil.SRA.w1.ox.ovary.inserts.xkxh.all.norm.bed.23-29.gz
		smmapper2=${smRNAINDIR}/Phil.SRA.${g}.ox.ovary.inserts.xkxh.all.norm.bed.23-29.gz #share the SRA norm.bed files
		demapper2=${INDIR}/Phil.DEG.${gt}.unox.ovary.PE.xkxh.${f}.mapper2.gz

	#total Ping-Pong
		[ ! -s ${OUT0}/${g}_${f}.total.pp8.out ] && echo "$script ${smmapper2} ${demapper2%.gz}.norm.bed.gz 2 ${OUTDIR}" >>${OUT0}/para.${g}_${f}.total.pp8.out 
		echo -e "`date` "+$ISO_8601"\ttotal Ping-Pong 8 analysis done..." >> $LOG
	done
done
fi
[ $? == 0 ] && \
touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8
STEP=$((STEP+1))


#generate master table for ppscore
masterOUT=${OUT0}/masterpp8score
[ ! -d ${masterOUT} ] && mkdir -p ${masterOUT}

touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master

if [ ! -f ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master ] 
then
	for f in "${FEATURE[@]}"
	do
		[ -f ${masterOUT}/SRA_norm.SRA_norm_${f}.nonnormalized.pp8score.txt ] && rm ${masterOUT}/SRA_norm.SRA_norm_${f}.nonnormalized.pp8score.txt
		[ -f ${masterOUT}/SRA_norm.DEG_${f}.nonnormalized.pp8score.txt ] && rm ${masterOUT}/SRA_norm.DEG_${f}.nonnormalized.pp8score.txt
		[ -f ${masterOUT}/DEG_${f}.SRA_norm.nonnormalized.pp8score.txt ] && rm ${masterOUT}/DEG_${f}.SRA_norm.nonnormalized.pp8score.txt
		[ -f ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.txt ] && rm ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.txt
		
		[ -f ${masterOUT}/SRA_norm.SRA_norm_${f}.normalized.pp8score.txt ] && rm ${masterOUT}/SRA_norm.SRA_norm_${f}.normalized.pp8score.txt
		[ -f ${masterOUT}/SRA_norm.DEG_${f}.normalized.pp8score.txt ] && rm ${masterOUT}/SRA_norm.DEG_${f}.normalized.pp8score.txt
		[ -f ${masterOUT}/DEG_${f}.SRA_norm.normalized.pp8score.txt ] && rm ${masterOUT}/DEG_${f}.SRA_norm.normalized.pp8score.txt
		[ -f ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.txt ] && rm ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.txt
		for g in "${GROUPGT[@]}"
		do
			OUTDIR=${OUT0}/${g}_${f}
			cut -f1,2 ${OUTDIR}/${g}_SRA_norm.${g}_SRA_norm.pp|awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}' >> ${masterOUT}/SRA_norm.SRA_norm_${f}.nonnormalized.pp8score.txt
			cut -f1,2 ${OUTDIR}/${g}_SRA_norm.${g}_DEG_${f}.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/SRA_norm.DEG_${f}.nonnormalized.pp8score.txt
			cut -f1,2 ${OUTDIR}/${g}_DEG_${f}.${g}_SRA_norm.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/DEG_${f}.SRA_norm.nonnormalized.pp8score.txt
			cut -f1,2 ${OUTDIR}/${g}_DEG_${f}.${g}_DEG_${f}.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.txt
	
			cut -f1,3 ${OUTDIR}/${g}_SRA_norm.${g}_SRA_norm.pp|awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}' >> ${masterOUT}/SRA_norm.SRA_norm_${f}.normalized.pp8score.txt
			cut -f1,3 ${OUTDIR}/${g}_SRA_norm.${g}_DEG_${f}.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/SRA_norm.DEG_${f}.normalized.pp8score.txt
			cut -f1,3 ${OUTDIR}/${g}_DEG_${f}.${g}_SRA_norm.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/DEG_${f}.SRA_norm.normalized.pp8score.txt
			cut -f1,3 ${OUTDIR}/${g}_DEG_${f}.${g}_DEG_${f}.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.txt
	
		done
		
	
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R_utils/R.source cast_master_table ${masterOUT}/SRA_norm.SRA_norm_${f}.nonnormalized.pp8score.txt ${masterOUT}/SRA_norm.SRA_norm_${f}.nonnormalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R_utils/R.source cast_master_table ${masterOUT}/SRA_norm.DEG_${f}.nonnormalized.pp8score.txt ${masterOUT}/SRA_norm.DEG_${f}.nonnormalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R_utils/R.source cast_master_table ${masterOUT}/DEG_${f}.SRA_norm.nonnormalized.pp8score.txt ${masterOUT}/DEG_${f}.SRA_norm.nonnormalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R_utils/R.source cast_master_table ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.txt ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.mastertable.txt
		
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R_utils/R.source cast_master_table ${masterOUT}/SRA_norm.SRA_norm_${f}.normalized.pp8score.txt ${masterOUT}/SRA_norm.SRA_norm_${f}.normalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R_utils/R.source cast_master_table ${masterOUT}/SRA_norm.DEG_${f}.normalized.pp8score.txt ${masterOUT}/SRA_norm.DEG_${f}.normalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R_utils/R.source cast_master_table ${masterOUT}/DEG_${f}.SRA_norm.normalized.pp8score.txt ${masterOUT}/DEG_${f}.SRA_norm.normalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R_utils/R.source cast_master_table ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.txt ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.mastertable.txt
	done
fi

[ $? == 0 ] && \
	touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master
STEP=$((STEP+1))