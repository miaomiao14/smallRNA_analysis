#! /bin/bash
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp6_T_ww_smRNA_vs_DEG.pl
smRNAINDIR=$1
degraINDIR=$2
OUT0=$3
g=$4
c=$5 #cpu
declare -a GROUPGT=("ago3MutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "aubMutsWW" "AubCDrescue" "AubWTrescue")

[ ! -d ${smRNAINDIR}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013 ] && mkdir -p ${smRNAINDIR}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013
#step 1
OUT=${OUT0}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013
for g in "${GROUPGT[@]}"
do
	[ ! -d ${smRNAINDIR}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013/${g} ] && mkdir -p ${smRNAINDIR}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013/${g}
	OUTDIR=${smRNAINDIR}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013/${g}
	LOG=${OUTDIR}/${g}.log
	smmapper2=${smRNAINDIR}/Phil.SRA.${g}.ox.ovary.inserts/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.gz
	demapper2=${OUTDIR}/Phil.DEG.${g}.ovary.PE.FLY_PIRNA_CLUSTER.xkxh.nta.mapper2.gz
	[ ! -f $demapper2 ] && \
		ln -s ${degraINDIR}/Phil.DEG.${g}.ovary.PE/bedIntersectWW/Phil.DEG.${g}.ovary.PE.x_rRNA.dm3.sorted.f0x40.noS.5p.all.bed.ntm.collapse.FLY_PIRNA_CLUSTER.nta.mapper2.gz $demapper2
	echo -e "`date` "+$ISO_8601"\tChanged the name of degradome transposon mapper2..." >> $LOG
	#length range: 23-29
	[ ! -f ${smmapper2%.gz}.23-29.gz ] && \
	${PIPELINE_DIRECTORY}/gzlenrangeselector.pl ${smmapper2} 23 29 >${smmapper2%.gz}.23-29 && \
	gzip ${smmapper2%.gz}.23-29
	echo -e "`date` "+$ISO_8601"\tSize select the smallRNA transposon mapper2 and gzip it..." >> $LOG
	#total Ping-Pong
	[ ! -s ${smRNAINDIR}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013/${g}.total.pp6.out ] && submitsge 8 $g ${OUT} "$script ${smmapper2%.gz}.23-29.gz ${demapper2} 2 ${smRNAINDIR}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013 >${smRNAINDIR}/pp6_FB_smRNAtrn_vs_degradomecluster_total_with_ppscore_11252013/${g}.total.pp6.out" 
	#echo -e "`date` "+$ISO_8601"\ttotal Ping-Pong analysis done..." >> $LOG
done