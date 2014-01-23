#! /bin/bash
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp8_ww_smRNA_vs_DEG.pl
smRNAINDIR=$1
degraINDIR=$2
OUT0=$3
declare -a FEATURE=("FLY_TRN_ALL_IN_CLUSTER")
#g=$5
#c=$6 #cpu
declare -a GROUPGT=("ago3MutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "aubMutsWW" "AubCDrescue" "AubWTrescue")

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
		smmapper2=${OUT0}/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.23-29.gz #share the SRA norm.bed files
		demapper2=${OUTDIR}/Phil.DEG.${g}.unox.ovary.PE.xkxh.${f}.mapper2.gz
		[ ! -f $demapper2 ] && \
			ln -s ${degraINDIR}/Phil.DEG.${g}.ovary.PE/bedIntersectWW/Phil.DEG.${g}.ovary.PE.x_rRNA.dm3.sorted.f0x40.noS.5p.all.bed.ntm.collapse.${f}.nta.mapper2.gz $demapper2
		[ ! -f $smmapper2 ] && \
			ln -s ${smRNAINDIR}/Phil.SRA.${g}.ox.ovary.inserts/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.23-29.gz ${smmapper2}	
		echo -e "`date` "+$ISO_8601"\tChanged the name of degradome transposon mapper2..." >> $LOG
			
	#length range: 23-29
	#[ ! -f ${smmapper2%.gz}.23-29.gz ] && \
	#${PIPELINE_DIRECTORY}/gzlenrangeselector.pl ${smmapper2} 23 29 >${smmapper2%.gz}.23-29 && \
	#gzip ${smmapper2%.gz}.23-29
		
	#convert mapper2 to norm.bed for pp8
		[ ! -s ${smmapper2%.gz}.23-29.norm.bed.gz ] && ${PIPELINE_DIRECTORY}/mapper2gznormbed.pl ${smmapper2} ${OUT0} && gzip ${smmapper2%.gz}.norm.bed
		[ ! -s ${demapper2%.gz}.norm.bed.gz ] && ${PIPELINE_DIRECTORY}/mapper2gznormbed.pl ${demapper2} ${OUTDIR} && gzip ${demapper2%.gz}.norm.bed
		
		echo -e "`date` "+$ISO_8601"\t convert the mapper2 format to norm.bed format done..." >> $LOG
		
		echo -e "`date` "+$ISO_8601"\tSize select the smallRNA transposon mapper2 and gzip it..." >> $LOG
	#total Ping-Pong
		[ ! -s ${OUT0}/${g}_${f}.total.pp8.out ] && submitsge 8 ${g}_${f} ${OUT0} "$script ${smmapper2%.gz}.norm.bed.gz ${demapper2%.gz}.norm.bed.gz 2 ${OUTDIR} >${OUT0}/${g}_${f}.total.pp8.out" 
		echo -e "`date` "+$ISO_8601"\ttotal Ping-Pong 8 analysis done..." >> $LOG
	done
done