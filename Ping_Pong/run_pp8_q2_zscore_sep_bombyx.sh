#!/bin/sh

#11/05/2013
#WEI WANG
#for Ping-Pong method paper, analyze the data from Yuki lab
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
INDIR=/home/wangw1/isilon_temp/BmN4
LOG=${INDIR}/LOG
indexFlag=$1

STEP=1
#run pp6 for total
touch .status.${STEP}.pp6forTotal
[ ! -f .status.${STEP}.pp6forTotal ] && \
totalBED=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.TOTAL.DMSO.ox.BmN4cell.inserts/Yuki.SRA.TOTAL.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.bed2 && \
submitsge 24 TOTALBEDPP6 ~/isilon_temp/BmN4/pp6_TOTAL "awk '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$4/\$5,\$6}' ${totalBED} >${totalBED%.bed2}.normalized.bed && /home/wangw1/git/smallRNA_analysis/pp6_ww_bed.pl -i ${totalBED%.bed2}.normalized.bed -o ~/isilon_temp/BmN4/Yuki.SRA.TOTAL.DMSO.ox.BmN4cell.inserts/pp6_total -f bedscore " && \
touch .status.${STEP}.pp6forTotal
STEP=$((STEP+1))

#run pp6 between Ago3 and Siwi
AGO3BED=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.bed2 
SIWIBED=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.bed2
touch .status.${STEP}.pp6forAgo3andSiwi
[ ! -f .status.${STEP}.pp6forAgo3andSiwi ] && \
submitsge 24 AGO3SIWIPP6 ~/isilon_temp/BmN4/pp6_TOTAL "awk '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$4/\$5,\$6}' ${AGO3BED} >${AGO3BED%.bed2}.normalized.bed && awk '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$4/\$5,\$6}' ${SIWIBED} >${SIWIBED%.bed2}.normalized.bed && \
/home/wangw1/git/smallRNA_analysis/pp6_ww_bed.pl -i ${AGO3BED%.bed2}.normalized.bed -j ${SIWIBED%.bed2}.normalized.bed -o ~/isilon_temp/BmN4/pp6_TOTAL -f bedscore" && \
touch .status.${STEP}.pp6forAgo3andSiwi
STEP=$((STEP+1))


#converting to mapper2 to norm.bed format,by mapper2normbed.pl
# run pp8_q2_ww1_zscore_sep_11052013.pl for UA_VA
declare -a GT=("Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell" "Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell" "Yuki.SRA.TOTAL.DMSO.ox.BmN4cell")
declare -a TARGETS=("GENE" "KNOWNTE" "ReASTE")
OUTDIR=/home/wangw1/isilon_temp/BmN4/pp8_q2_repeat

if [ ! -f .status.${STEP}.pp8_UA_VA_repeat ] 
then
for t in ${TARGETS[@]}
do
	fasta=/home/wangw1/pipeline_bm/common/silkgenome.fa
	A=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.${t}.AS.norm.bed.gz
	B=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.${t}.S.norm.bed.gz
	jobname=${t}_Ago3AS_SiwiS.pp8.q2
	OUT=${OUTDIR}/${t}_Ago3AS_SiwiS
	[ ! -d ${OUT} ] && mkdir -p ${OUT}
	/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "/home/wangw1/git/smallRNA_analysis/Ping_Pong/pp8_q2_ww1_zscore_sep_11052013.pl ${A} ${B} 2 bombyx ${OUT} ${indexFlag} >${OUTDIR}/${t}_Ago3AS_SiwiS.pp8.q2.UA_VA.log"

	A=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.${t}.S.norm.bed.gz
	B=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.${t}.AS.norm.bed.gz
	jobname=${t}_Ago3S_SiwiAS.pp8.q2
	OUT=${OUTDIR}/${t}_Ago3S_SiwiAS
	[ ! -d ${OUT} ] && mkdir -p ${OUT}
	/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "/home/wangw1/git/smallRNA_analysis/Ping_Pong/pp8_q2_ww1_zscore_sep_11052013.pl ${A} ${B} 2 bombyx ${OUT} ${indexFlag} >${OUTDIR}/${t}_Ago3S_SiwiAS.pp8.q2.UA_VA.log"	
done
fi
[ $? == 0 ] && \
touch .status.${STEP}.pp8_UA_VA_repeat
STEP=$((STEP+1))

declare -a PPPAIR=("Ago3AS_SiwiS" "Ago3S_SiwiAS")
declare -a TARGETS=("KNOWNTE" "ReASTE")
OUTDIR=/home/wangw1/isilon_temp/BmN4/pp8_q2_repeat_summary
[ ! -d ${OUTDIR} ] && mkdir $OUTDIR

if [ ! -f .status.${STEP}.pp8_UA_VA_repeat_summary ] 
then
for t in ${TARGETS[@]}
do
	for pp in ${PPPAIR[@]}
	do
		fn=/home/wangw1/isilon_temp/BmN4/pp8_q2_repeat/${t}_${pp}
		file=/home/wangw1/isilon_temp/BmN4/pp8_q2_repeat_summary/${t}_${pp}
		[ -f ${file}.pair.count.txt ] &&  rm ${file}.pair.count.txt
		for i in `ls ${fn}/*.UA_VA.zscore.out`
		do
		[ -f ${i} ] && /home/wangw1/git/smallRNA_analysis/Ping_Pong/UA_VA_rawscore.pl ${i} $OUTDIR >> ${file}.pair.count.txt
		done
		#sort -k1,1 -k2,2 -k3,3 -k4,4 $file.pair.count.txt | uniq >${file}.pair.count.uniq.txt
		#RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va ${file}.pair.count.uniq.txt ${t}_${pp} ${OUTDIR}
		[ -f ${file}.pair.count.uniq.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va_color ${file}.pair.count.uniq.txt ${t}_${pp} ${OUTDIR}
	done
	
	
done
fi
[ $? == 0 ] && \
touch .status.${STEP}.pp8_UA_VA_repeat_summary
STEP=$((STEP+1))


declare -a PPPAIR=("Ago3AS_SiwiS" "Ago3S_SiwiAS")
declare -a TARGETS=("KNOWNTE" "ReASTE")
OUTDIR=/home/wangw1/isilon_temp/BmN4/pp8_q2_repeat_summary_from_PP
[ ! -d ${OUTDIR} ] && mkdir $OUTDIR

if [ ! -f .status.${STEP}.pp8_UA_VA_repeat_summary_from_PP ] 
then
for t in ${TARGETS[@]}
do
	for pp in ${PPPAIR[@]}
	do
		fn=/home/wangw1/isilon_temp/BmN4/pp8_q2_repeat/${t}_${pp}
		file=/home/wangw1/isilon_temp/BmN4/pp8_q2_repeat_summary_from_PP/${t}_${pp}
		[ -f ${file}.pair.count.txt ] &&  rm ${file}.pair.count.txt
		for i in `ls ${fn}/*.VA.pp`
		do
			filename=${i##*/}
			pairname=`basename ${filename} .VA.pp`
			[ -f $i ] && awk -v gt=${pairname} '{OFS="\t"}{print gt,$0}' ${i} >${i}.gt && \
			/home/wangw1/git/smallRNA_analysis/Ping_Pong/UA_VA_rawscore_from_ppscore.pl ${i}.gt $OUTDIR >> ${file}.pair.count.txt 
		
		done
		#sort -k1,1 -k2,2 -k3,3 -k4,4 $file.pair.count.txt | uniq >${file}.pair.count.uniq.txt
		#RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va ${file}.pair.count.uniq.txt ${t}_${pp} ${OUTDIR}
		[ -f ${file}.pair.count.uniq.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va_from_ppscore_color ${file}.pair.count.uniq.txt ${t}_${pp} ${OUTDIR}
	done
	
	
done
fi
[ $? == 0 ] && \
touch .status.${STEP}.pp8_UA_VA_repeat_summary
STEP=$((STEP+1))



echo -e "`date` "+$ISO_8601"\tDraw phasing analysis..." >> $LOG

OUTDIR=/home/wangw1/isilon_temp/BmN4/phasing
[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}
touch .status.${STEP}.transposon_piRNA.phasing
if [ ! -f .status.${STEP}.transposon_piRNA.phasing ] 
then
#paraFile=${OUTDIR13}/${RANDOM}.piRNAphasing.para

for i in `ls ${INDIR}/*.TOTAL.*.inserts/*.norm.bed.gz`
do
	inputfile=${i##*/}
	samplenamepart=${inputfile#Yuki.SRA.*}
	samplename=${samplenamepart%*.norm.bed.gz}
	sample=${samplename/DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh./}
	/home/wangw1/bin/submitsge 8 ${sample} $OUTDIR "${PIPELINE_DIRECTORY}/run_distance_analysis.sh -i ${i} -o $OUTDIR -t normbed" 
done
fi
[ $? == 0 ] && \
touch .status.${STEP}.transposon_piRNA.phasing

echo -e "`date` "+$ISO_8601"\tDraw phasing analysis done" >> $LOG
#zip the mapper2 format, run pp6
#Ago3IP S: SiwiIP AS, Ago3IP AS: SiwiIP S, Ago3IP total: SiwiIP total	