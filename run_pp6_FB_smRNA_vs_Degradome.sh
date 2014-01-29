#! /bin/bash
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp6_T_ww_smRNA_vs_DEG.pl
smRNAINDIR=$1
degraINDIR=$2
OUT0=$3
c=$4 #cpu
g=$5 
declare -a GROUPGT=("ago3MutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "aubMutsWW" "AubCDrescue" "AubWTrescue")
declare -a FEATURE=("FLY_TRN_ALL" "FLY_TRN_ALL_IN_CLUSTER" "FLY_TRN_ALL_OUT_CLUSTER")


[ ! -d ${OUT0}/pp6_FB_smRNA_vs_degradome_FB_withppscore_11252013 ] && mkdir -p ${OUT0}/pp6_FB_smRNA_vs_degradome_FB_withppscore_11252013
#step 1
#for g in "${GROUPGT[@]}"
#do
for g in "${GROUPGT[@]}"
do
	[ ! -d ${OUT0}/pp6_FB_smRNA_vs_degradome_FB_withppscore_11252013/${g} ] && mkdir -p ${OUT0}/pp6_FB_smRNA_vs_degradome_FB_withppscore_11252013/${g}
	OUTDIR=${OUT0}/pp6_FB_smRNA_vs_degradome_FB_withppscore_11252013/${g}
	LOG=${OUTDIR}/${g}.log
	smmapper2=${smRNAINDIR}/Phil.SRA.${g}.ox.ovary.inserts/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.gz
	demapper2=${OUT0}/pp6_FB_smRNA_vs_degradome_FB_withppscore_11252013/Phil.DEG.${g}.ovary.PE.FLY_TRANSPOSON_ALL.xkxh.nta.mapper2.gz
	[ ! -f $demapper2 ] && \
		ln -s ${degraINDIR}/Phil.DEG.${g}.ovary.PE/bedIntersectWW/Phil.DEG.${g}.ovary.PE.x_rRNA.dm3.sorted.f0x40.noS.5p.all.bed.ntm.collapse.FLY_TRANSPOSON_ALL.nta.mapper2.gz $demapper2
	echo -e "`date` "+$ISO_8601"\tChanged the name of degradome transposon mapper2..." >> $LOG
	#length range: 23-29
	[ ! -f ${smmapper2%.gz}.23-29.gz ] && \
	${PIPELINE_DIRECTORY}/gzlenrangeselector.pl ${smmapper2} 23 29 >${smmapper2%.gz}.23-29 && \
	gzip ${smmapper2%.gz}.23-29
	echo -e "`date` "+$ISO_8601"\tSize select the smallRNA transposon mapper2 and gzip it..." >> $LOG
	#total Ping-Pong
	#[ ! -s ${smRNAINDIR}/pp6_FB_smRNA_vs_degradome/${g}.total.pp6.out ] && $script ${smmapper2%.gz}.23-29.gz ${demapper2} 2 ${smRNAINDIR}/pp6_FB_smRNA_vs_degradome >${smRNAINDIR}/pp6_FB_smRNA_vs_degradome/${g}.total.pp6.out && \
	#echo -e "`date` "+$ISO_8601"\ttotal Ping-Pong analysis done..." >> $LOG
	
	[ ! -s ${OUTDIR}/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.23-29.roo ] && ${PIPELINE_DIRECTORY}/FB.pl ${smmapper2%.gz}.23-29.gz ${OUTDIR}
	[ ! -s ${OUTDIR}/Phil.DEG.${g}.ovary.PE.FLY_TRANSPOSON_ALL.xkxh.nta.mapper2.roo ] && ${PIPELINE_DIRECTORY}/FB.pl ${demapper2} ${OUTDIR} && \
	echo -e "`date` "+$ISO_8601"\tsplit transposon mapper2 to trns..." >> $LOG
	paraFile=${OUTDIR}/${g}.pp6.para
	for j in `ls -1 ${OUTDIR}/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.23-29.*`
	do
	T1=${j##*mapper2.23-29.}
	[ -s $OUTDIR/${g}.FB.pp6.temp ] && rm $OUTDIR/${g}.FB.pp6.temp
	[ ! -d $OUTDIR/${T1} ] && mkdir -p $OUTDIR/${T1} 
	echo -ne  " T1=${j##*mapper2.23-29.} && " >> ${paraFile}
	echo -ne " $script ${j}  ${OUTDIR}/Phil.DEG.${g}.ovary.PE.FLY_TRANSPOSON_ALL.xkxh.nta.mapper2.${T1} 2 $OUTDIR/${T1} $T1 >> $OUTDIR/${T1}/${g}.FB.${T1}.pp6.out && " >> ${paraFile}
	echo -e " $script ${j} ${OUTDIR}/Phil.DEG.${g}.ovary.PE.FLY_TRANSPOSON_ALL.xkxh.nta.mapper2.${T1} 2 $OUTDIR/${T1} $T1 >> $OUTDIR/${g}.FB.pp6.temp  " >> ${paraFile}
	done
	if [[ ! -f ${paraFile}.completed ]] || [[ -f $paraFile.failed_commands ]]
	then
	
		ParaFly -c $paraFile -CPU $c -failed_cmds $paraFile.failed_commands
	fi
	echo -e "`date` "+$ISO_8601"\ttrn PP done..." >> $LOG
	awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${OUTDIR}/${g}.FB.pp6.temp > ${OUTDIR}/${g}.FB.pp6
done