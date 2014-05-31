#!/bin/bash -x
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp8_ww_smRNA_vs_DEG.pl
INDIR=/home/wangw1/data/projects/cd/smRNA/pp8_SRADEG_datasets
smRNAINDIR=${INDIR}/smRNA/jia_pipeline_results
degraINDIR=${INDIR}/degradome/pipeline_output_12262013
OUT0=/home/wangw1/data/projects/cd/smRNA/pp8_smRNAall_vs_degradometrnoutcluster_total_05312014

#note: use window size 16 for complementarity 16 (updated into script)
declare -a FEATURE=("FLY_TRN_ALL" "FLY_TRN_ALL_IN_CLUSTER" "FLY_TRN_ALL_OUT_CLUSTER")
#g=$5
#c=$6 #cpu
declare -a GROUPGT=("ago3MutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "aubMutsWW" "AubCDrescue" "AubWTrescue")
#declare -a GROUPGT=("ago3MutsWW")

STEP=1

[ ! -d ${OUT0} ] && mkdir -p ${OUT0}
#step 1
#OUT=${OUT0}
#touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8
if [ ! -f ${OUT0}/.status.${STEP}.SRA_DEG.pp8 ] 
then
for g in "${GROUPGT[@]}"
do
	
	for f in "${FEATURE[@]}"
	do
		OUTDIR=${OUT0}/${g}_${f}
		[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
		
		LOG=${OUTDIR}/${g}.log	
		smmapper2=${INDIR}/Phil.SRA.${g}.ox.ovary.inserts.xkxh.norm.bed.gz #share the SRA norm.bed files
		demapper2=${INDIR}/Phil.DEG.${g}.unox.ovary.PE.xkxh.${f}.mapper2.gz

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
		[ -f ${masterOUT}/SRA_transposon.SRA_transposon_${f}.nonnormalized.pp8score.txt ] && rm ${masterOUT}/SRA_transposon.SRA_transposon_${f}.nonnormalized.pp8score.txt
		[ -f ${masterOUT}/SRA_transposon.DEG_${f}.nonnormalized.pp8score.txt ] && rm ${masterOUT}/SRA_transposon.DEG_${f}.nonnormalized.pp8score.txt
		[ -f ${masterOUT}/DEG_${f}.SRA_transposon.nonnormalized.pp8score.txt ] && rm ${masterOUT}/DEG_${f}.SRA_transposon.nonnormalized.pp8score.txt
		[ -f ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.txt ] && rm ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.txt
		
		[ -f ${masterOUT}/SRA_transposon.SRA_transposon_${f}.normalized.pp8score.txt ] && rm ${masterOUT}/SRA_transposon.SRA_transposon_${f}.normalized.pp8score.txt
		[ -f ${masterOUT}/SRA_transposon.DEG_${f}.normalized.pp8score.txt ] && rm ${masterOUT}/SRA_transposon.DEG_${f}.normalized.pp8score.txt
		[ -f ${masterOUT}/DEG_${f}.SRA_transposon.normalized.pp8score.txt ] && rm ${masterOUT}/DEG_${f}.SRA_transposon.normalized.pp8score.txt
		[ -f ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.txt ] && rm ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.txt
		for g in "${GROUPGT[@]}"
		do
			OUTDIR=${OUT0}/${g}_${f}
			cut -f1,2 ${OUTDIR}/${g}_SRA_transposon.${g}_SRA_transposon.pp|awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}' >> ${masterOUT}/SRA_transposon.SRA_transposon_${f}.nonnormalized.pp8score.txt
			cut -f1,2 ${OUTDIR}/${g}_SRA_transposon.${g}_DEG_${f}.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/SRA_transposon.DEG_${f}.nonnormalized.pp8score.txt
			cut -f1,2 ${OUTDIR}/${g}_DEG_${f}.${g}_SRA_transposon.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/DEG_${f}.SRA_transposon.nonnormalized.pp8score.txt
			cut -f1,2 ${OUTDIR}/${g}_DEG_${f}.${g}_DEG_${f}.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.txt
	
			cut -f1,3 ${OUTDIR}/${g}_SRA_transposon.${g}_SRA_transposon.pp|awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}' >> ${masterOUT}/SRA_transposon.SRA_transposon_${f}.normalized.pp8score.txt
			cut -f1,3 ${OUTDIR}/${g}_SRA_transposon.${g}_DEG_${f}.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/SRA_transposon.DEG_${f}.normalized.pp8score.txt
			cut -f1,3 ${OUTDIR}/${g}_DEG_${f}.${g}_SRA_transposon.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/DEG_${f}.SRA_transposon.normalized.pp8score.txt
			cut -f1,3 ${OUTDIR}/${g}_DEG_${f}.${g}_DEG_${f}.pp |awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}'>> ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.txt
	
		done
		
	
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/SRA_transposon.SRA_transposon_${f}.nonnormalized.pp8score.txt ${masterOUT}/SRA_transposon.SRA_transposon_${f}.nonnormalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/SRA_transposon.DEG_${f}.nonnormalized.pp8score.txt ${masterOUT}/SRA_transposon.DEG_${f}.nonnormalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/DEG_${f}.SRA_transposon.nonnormalized.pp8score.txt ${masterOUT}/DEG_${f}.SRA_transposon.nonnormalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.txt ${masterOUT}/DEG_${f}.DEG_${f}.nonnormalized.pp8score.mastertable.txt
		
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/SRA_transposon.SRA_transposon_${f}.normalized.pp8score.txt ${masterOUT}/SRA_transposon.SRA_transposon_${f}.normalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/SRA_transposon.DEG_${f}.normalized.pp8score.txt ${masterOUT}/SRA_transposon.DEG_${f}.normalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/DEG_${f}.SRA_transposon.normalized.pp8score.txt ${masterOUT}/DEG_${f}.SRA_transposon.normalized.pp8score.mastertable.txt
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.txt ${masterOUT}/DEG_${f}.DEG_${f}.normalized.pp8score.mastertable.txt
	done
fi

[ $? == 0 ] && \
	touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master
STEP=$((STEP+1))