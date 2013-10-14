#! /bin/bash
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp6_T_ww_smRNA_vs_DEG.pl
smRNAINDIR=$1
degraINDIR=$2
g=$3
c=$4 #cpu
declare -a GROUPGT=("ago3MutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "aubMutsWW" "AubCDrescue" "AubWTrescue")

[ ! -d ${smRNAINDIR}/pp6_FB_smRNAuniq_vs_degradomeuniq_total_with_ppscore ] && mkdir -p ${smRNAINDIR}/pp6_FB_smRNAuniq_vs_degradomeuniq_total_with_ppscore
#step 1
for g in "${GROUPGT[@]}"
do
	[ ! -d ${smRNAINDIR}/pp6_FB_smRNAuniq_vs_degradomeuniq_total_with_ppscore/${g} ] && mkdir -p ${smRNAINDIR}/pp6_FB_smRNAuniq_vs_degradomeuniq_total_with_ppscore/${g}
	OUTDIR=${smRNAINDIR}/pp6_FB_smRNAuniq_vs_degradomeuniq_total_with_ppscore/${g}
	LOG=${OUTDIR}/${g}.log
	smmapper2=${smRNAINDIR}/Phil.SRA.${g}.ox.ovary.inserts/Phil.SRA.${g}.ox.ovary.inserts.uniqmap.xkxh.transposon.mapper2.gz
	demapper2=${degraINDIR}/Phil.DEG.${g}.ovary/bedIntersectWW/Phil.DEG.${g}.ovary.x_rRNA.dm3.sorted.f0x40.noS.5p.unique.bed.ntm.collapse.xkxh.FLY_TRANSPOSON_ALL.nta.mapper2.gz
	[ ! -f $demapper2 ] && \
		mv ${degraINDIR}/Phil.DEG.${g}.ovary/bedIntersectWW/Phil.DEG.${g}.ovary.x_rRNA.dm3.sorted.f0x40.noS.5p.unique.bed.ntm.collapse.FLY_TRANSPOSON_ALL.nta.mapper2.gz $demapper2
	echo -e "`date` "+$ISO_8601"\tChanged the name of degradome transposon mapper2..." >> $LOG
	#length range: 23-29
	[ ! -f ${smmapper2%.gz}.23-29.gz ] && \
	${PIPELINE_DIRECTORY}/gzlenrangeselector.pl ${smmapper2} 23 29 >${smmapper2%.gz}.23-29 && \
	gzip ${smmapper2%.gz}.23-29
	echo -e "`date` "+$ISO_8601"\tSize select the smallRNA transposon mapper2 and gzip it..." >> $LOG
	#total Ping-Pong
	[ ! -s ${smRNAINDIR}/pp6_FB_smRNAuniq_vs_degradomeuniq_total_with_ppscore/${g}.total.pp6.out ] && submitsge 8 $g "$script ${smmapper2%.gz}.23-29.gz ${demapper2} 2 ${smRNAINDIR}/pp6_FB_smRNAuniq_vs_degradomeuniq_total_with_ppscore >${smRNAINDIR}/pp6_FB_smRNAuniq_vs_degradomeuniq_total_with_ppscore/${g}.uniqmap.total.pp6.out" 
	#echo -e "`date` "+$ISO_8601"\ttotal Ping-Pong analysis done..." >> $LOG
done