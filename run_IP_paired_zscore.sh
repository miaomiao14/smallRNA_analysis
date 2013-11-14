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

INDIR=$1 #this is the folder store all pipeline results outmost folders 
#/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log

STEP=1

declare -a GROUPGT=("AubIP_ago3cdwt_ox" "AubIP_ago3cdwt_unox" \
"AubIP_ago3cdw1_ox" "AubIP_ago3cdw1_unox" \
)

declare -a AubIP_ago3cdwt_ox=("Phil.AubIP.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.AubIP.aubvasAgo3WTrescue.ox.ovary.inserts")
declare -a AubIP_ago3cdwt_unox=("Phil.AubIP.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.AubIP.aubvasAgo3WTrescue.unox.ovary.inserts")

declare -a AubIP_ago3cdw1_ox=("Phil.AubIP.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.AubIP.w1.ox.ovary.inserts")
declare -a AubIP_ago3cdw1_unox=("Phil.AubIP.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.AubIP.w1.unox.ovary.inserts")

echo -e "`date` "+$ISO_8601"\tDraw paired transposon piRNA zscore scatterplot..." >> $LOG
OUTDIR10=${INDIR}/transposon_piRNA/paired_zscore_scatterplot
[ ! -d $OUTDIR10 ] && mkdir -p ${OUTDIR10}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.pairedZscore ] && \
paraFile=${OUTDIR10}/${RANDOM}.pairedZscore.para && \

for g in "${GROUPGT[@]}"
do
	SUBGROUP="$g[@]"
	[ -f ${OUTDIR10}/${g}.FB.zscore ] && rm ${OUTDIR10}/${g}.FB.zscore
	[ -f ${OUTDIR10}/${g}.total.zscore ] && rm ${OUTDIR10}/${g}.total.zscore
	count=1
	for t in ${!SUBGROUP}
	do
		cat ${INDIR}/pp6_FB/${t}/${t}.FB.pp6.out| awk -v gt=$t -v rank=$count '{OFS="\t"}{print gt,$1,$2,$3,$4,rank}' >> ${OUTDIR10}/${g}.FB.zscore
		count=$(($count+1))
		totalZscore=`cat ${INDIR}/pp6_FB/${t}/${t}.total.pp6.out|cut -f2`
		TZ=`printf "%0.2f" $totalZscore`
		echo -e ${t}"\t"${TZ} >>${OUTDIR10}/${g}.total.zscore
	done
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_zscore_FB_scatterplot ${OUTDIR10}/${g}.FB.zscore ${OUTDIR10}/${g}.total.zscore $g $OUTDIR10 " >>${paraFile}
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.transposon_piRNA.pairedZscore
STEP=$((STEP+1))



#run UA_VA analysis
#separate sense and antisense from normbed
#Phil.AubIP.AubCDrescue.unox.ovary.inserts.xkxh.norm.bed.transposons.gz
echo -e "`date` "+$ISO_8601"\t1U10A and 1V10A analysis..." >> $LOG
echo -e "`date` "+$ISO_8601"\tseparate the sense and antisense transposon mappers in norm.bed format..." >> $LOG
#touch ${OUT}/.status.transposon_piRNA.S_AS
OUTDIR=${INDIR}/transposon_piRNA/UA_VA
[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}
[ ! -f ${OUT}/.status.transposon_piRNA.S_AS ] && \
paraFile=${OUTDIR}/UA_VA.${RANDOM}.para && \
for i in `find ${INDIR} -name "*.inserts" -maxdepth 1`
do
	name=${i##*/}
	echo -ne "gunzip -f ${i}/${name}.xkxh.transposon.mapper2.gz && mapper2normbed.pl ${i}/${name}.xkxh.transposon.mapper2 >${i}/${name}.xkxh.norm.bed.transposons && ">>$paraFile;
	echo -ne "grep -w sense ${i}/${name}.xkxh.norm.bed.transposons >${i}/${name}.xkxh.norm.bed.transposons.sense && ">>$paraFile;
	
	echo -ne "gzip -f ${i}/${name}.xkxh.norm.bed.transposons.sense && ">>$paraFile;
	echo -ne "grep -w antisense ${i}/${name}.xkxh.norm.bed.transposons >${i}/${name}.xkxh.norm.bed.transposons.antisense && ">>$paraFile;
	echo -ne "gzip -f ${i}/${name}.xkxh.norm.bed.transposons.antisense && ">>$paraFile;
	echo -e "gzip -f ${i}/${name}.xkxh.transposon.mapper2 && rm ${i}/${name}.xkxh.norm.bed.transposons ">>$paraFile;
done
[ $? == 0 ] && \
	/home/wangw1/bin/submitsge 24 paraflyrun ${OUTDIR} "ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands"
[ $? == 0 ] && \
	touch ${OUT}/.status.transposon_piRNA.S_AS 
echo -e "`date` "+$ISO_8601"\tseparate the sense and antisense transposon mappers in norm.bed format done!" >> $LOG	


echo -e "`date` "+$ISO_8601"\trun pp8 UA and VA..." >> $LOG
declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "ago3Hets" "aubHets" "qinHets" "nosAgo3CDrescue" "nosAgo3WTrescue")
declare -a OX=("ox" "unox")
declare -a UNIQ=("uniq" "shared")

#uniq
#shared
indexFlag=1 #to indicate we need to build the index or not
	

[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.UA_VA ] && \
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do
		for s in ${UNIQ[@]}
		do
			A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.sense.gz
			B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.antisense.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}${s}_${o}_Ago3S_AubAS.pp8.q2 && \
			jOUT=${OUTDIR}/${t}${s}_${o}_Ago3S_AubAS && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "/home/wangw1/git/smallRNA_analysis/Ping_Pong/pp8_q2_ww1_zscore_sep_11052013.pl ${A} ${B} 2 fly ${jOUT} ${indexFlag} >${jOUT}/${t}${s}_${o}_Ago3S_AubAS.pp8.q2.UA_VA.log"
#			
			A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.antisense.gz
			B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.sense.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}${s}_${o}_Ago3AS_AubS.pp8.q2 && \
			jOUT=${OUTDIR}/${t}${s}_${o}_Ago3AS_AubS && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "/home/wangw1/git/smallRNA_analysis/Ping_Pong/pp8_q2_ww1_zscore_sep_11052013.pl ${A} ${B} 2 fly ${jOUT} ${indexFlag} >${jOUT}/${t}${s}_${o}_Ago3AS_AubS.pp8.q2.UA_VA.log"			
		done
	done
done
#total
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do

			A=${INDIR}/Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.sense.gz
			B=${INDIR}/Phil.SRA.AubIP.${t}.${o}.ovary.inserts/Phil.SRA.AubIP.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.antisense.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}_${o}_Ago3S_AubAS.pp8.q2 && \
			jOUT=${OUTDIR}/${t}_${o}_Ago3S_AubAS && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "/home/wangw1/git/smallRNA_analysis/Ping_Pong/pp8_q2_ww1_zscore_sep_11052013.pl ${A} ${B} 2 fly ${jOUT} ${indexFlag} >${jOUT}/${t}_${o}_Ago3S_AubAS.pp8.q2.UA_VA.log"
			
			A=${INDIR}/Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.antisense.gz
			B=${INDIR}/Phil.SRA.AubIP.${t}.${o}.ovary.inserts/Phil.SRA.AubIP.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.sense.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}_${o}_Ago3AS_AubS.pp8.q2 && \
			jOUT=${OUTDIR}/${t}_${o}_Ago3AS_AubS && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "/home/wangw1/git/smallRNA_analysis/Ping_Pong/pp8_q2_ww1_zscore_sep_11052013.pl ${A} ${B} 2 fly ${jOUT} ${indexFlag} >${jOUT}/${t}_${o}_Ago3AS_AubS.pp8.q2.UA_VA.log"			

	done
done
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.transposon_piRNA.UA_VA


#IPed piRNA abundance normalzation



