#! /bin/bash
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp6_T_ww_smRNA_vs_DEG.pl
smRNAINDIR=$1
degraINDIR=$2
g=$3
c=$4 #cpu
declare -a GROUPGT=("ago3MutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "aubMutsWW" "AubCDrescue" "AubWTrescue")

[ ! -d ${smRNAINDIR}/pp6_FB_smRNA_vs_degradome ] && mkdir -p ${smRNAINDIR}/pp6_FB_smRNA_vs_degradome
#step 1
#for g in "${GROUPGT[@]}"
#do
	[ ! -d ${smRNAINDIR}/pp6_FB_smRNA_vs_degradome/${g} ] && mkdir -p ${smRNAINDIR}/pp6_FB_smRNA_vs_degradome/${g}
	OUTDIR=${smRNAINDIR}/pp6_FB_smRNA_vs_degradome/${g}
	smmapper2=${smRNAINDIR}/Phil.SRA.${g}.ox.ovary.inserts/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.gz
	demapper2=${degraINDIR}/Phil.DEG.${g}.ovary.PE/bedIntersectWW/Phil.DEG.${g}.ovary.x_rRNA.dm3.sorted.f0x40.noS.5p.unique.bed.ntm.collapse.FLY_TRANSPOSON_ALL.nta.mapper2.gz
	
	#length range: 23-29
	${PIPELINE_DIRECTORY}/gzlenrangeselector.pl ${smmapper2} 23 29 >${smmapper2%.gz}.23-29
	gzip ${smmapper2%.gz}.23-29
	#total Ping-Pong
	[ ! -f ${smRNAINDIR}/pp6_FB_smRNA_vs_degradome/${g}.total.pp6.out ] && $script ${smmapper2%.gz}.23-29.gz ${demapper2} 2 ${smRNAINDIR}/pp6_FB_smRNA_vs_degradome >${smRNAINDIR}/pp6_FB_smRNA_vs_degradome/${g}.total.pp6.out &&
	
	
	${PIPELINE_DIRECTORY}/FB.pl ${smmapper2%.gz}.23-29.gz ${OUTDIR} &&
	${PIPELINE_DIRECTORY}/FB.pl ${demapper2} ${OUTDIR}	&& 
	paraFile=${OUTDIR}/${g}.pp6.para
	for j in `ls -1 ${OUTDIR}/Phil.SRA.${g}.ox.ovary.inserts.xkxh.transposon.mapper2.23-29.*`
	do
	T1=${j##*mapper2.}
	echo -ne  " T1=\${j##*mapper2.} && \" >> ${paraFile}
	echo -ne " $script ${j}  ${OUTDIR}/Phil.DEG.${g}.ovary.x_rRNA.dm3.sorted.f0x40.noS.5p.unique.bed.ntm.collapse.FLY_TRANSPOSON_ALL.nta.mapper2.${T1} 2 ${OUTDIR} $T1 >> $OUTDIR/${g}.FB.${T1}.pp6.out && \" >> ${paraFile}
	echo -e " $script ${j} ${OUTDIR}/Phil.DEG.${g}.ovary.x_rRNA.dm3.sorted.f0x40.noS.5p.unique.bed.ntm.collapse.FLY_TRANSPOSON_ALL.nta.mapper2.${T1} 2 ${OUTDIR} $T1 >> $OUTDIR/${g}.FB.pp6.temp  " >> ${paraFile}
	done
	if [[ ! -f \${paraFile}.completed ]] || [[ -f \$paraFile.failed_commands ]]
	then
	
		ParaFly -c \$paraFile -CPU $c -failed_cmds \$paraFile.failed_commands
	fi
	awk '{OFS="\t"}{print $1,$2,$3,$4}' ${OUTDIR}/${g}.FB.pp6.temp > ${OUTDIR}/${g}.FB.pp6
#done