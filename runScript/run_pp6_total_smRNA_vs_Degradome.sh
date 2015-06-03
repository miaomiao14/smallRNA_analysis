#! /bin/bash
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp6_T_ww_smRNA_vs_DEG.pl
smRNAINDIR=/home/wangw1/isilon_temp/smRNA/jia_pipeline_results
degraINDIR=/home/wangw1/isilon_temp/degradome/pipeline_output_12262013
OUT0=/home/wangw1/isilon_temp/smRNA/pp6_smRNAtrn_vs_degradometrnoutcluster_total_01292014
#FEATURE=$4
#g=$5
#c=$6 #cpu
declare -a GROUPGT=("ago3MutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "aubMutsWW" "AubCDrescue" "AubWTrescue")
declare -a FEATURE=("FLY_TRN_ALL" "FLY_TRN_ALL_IN_CLUSTER" "FLY_TRN_ALL_OUT_CLUSTER")
[ ! -d ${OUT0} ] && mkdir -p ${OUT0}
#step 1
#OUT=${OUT0}
for g in "${GROUPGT[@]}"
do
	for f in "${FEATURE[@]}"
	do
	OUTDIR=${OUT0}/${g}_${f}
	[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
	
	LOG=${OUTDIR}/${g}.log
	smmapper2=${smRNAINDIR}/Phil.SRA.${g}.ox.ovary.inserts/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.gz
	demapper2=${OUTDIR}/Phil.DEG.${g}.ovary.PE.${f}.xkxh.nta.mapper2.gz
	[ ! -f $demapper2 ] && \
		ln -s ${degraINDIR}/Phil.DEG.${g}.ovary.PE/bedIntersectWW/Phil.DEG.${g}.ovary.PE.x_rRNA.dm3.sorted.f0x40.noS.5p.all.bed.ntm.collapse.${f}.nta.mapper2.gz $demapper2
	echo -e "`date` "+$ISO_8601"\tChanged the name of degradome transposon mapper2..." >> $LOG
	#length range: 23-29
	[ ! -f ${smmapper2%.gz}.23-29.gz ] && \
	${PIPELINE_DIRECTORY}/gzlenrangeselector.pl ${smmapper2} 23 29 >${smmapper2%.gz}.23-29 && \
	gzip ${smmapper2%.gz}.23-29
	echo -e "`date` "+$ISO_8601"\tSize select the smallRNA transposon mapper2 and gzip it..." >> $LOG
	#total Ping-Pong
	[ ! -s ${OUT0}/${g}.total.pp6.out ] && submitsge 8 $g ${OUT0} "$script ${smmapper2%.gz}.23-29.gz ${demapper2} 2 ${OUTDIR} >${OUT0}/${g}.total.pp6.out" 
	#echo -e "`date` "+$ISO_8601"\ttotal Ping-Pong analysis done..." >> $LOG
	done
done