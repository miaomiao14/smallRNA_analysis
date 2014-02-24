#!/bin/sh

#11/05/2013
#WEI WANG
#for Ping-Pong method paper, analyze the data from Yuki lab
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
INDIR=/home/wangw1/isilon_temp/BmN4
LOG=${INDIR}/LOG
indexFlag=$1
script=${PIPELINE_DIRECTORY}/g1Ut10A/pp8_q2_ww1_zscore_sep_02102014.pl
STEP=1

#converting to mapper2 to norm.bed format,by mapper2normbed.pl
# run pp8_q2_ww1_zscore_sep_11052013.pl for UA_VA
declare -a GT=("Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell" "Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell")
declare -a TARGETS=("KNOWNTE")
OUTDIR=/home/wangw1/isilon_temp/BmN4/pp8_UA_VA_cistrans
[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
if [ ! -f .status.${STEP}.pp8_UA_VA_cistrans ] 
then
for t in ${TARGETS[@]}
do
	fasta=/home/wangw1/pipeline_bm/common/silkgenome.formatted.fa
	A=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.${t}.AS.norm.bed.gz
	B=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.${t}.S.norm.bed.gz
	jobname=${t}_Ago3AS_SiwiS.pp8.q2
	OUT=${OUTDIR}/${t}_Ago3AS_SiwiS
	[ ! -d ${OUT} ] && mkdir -p ${OUT}
	/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "${script} -i ${A} -j ${B} -n 2 -s bombyx -o ${OUT} -d ${indexFlag} -f normbed >${OUTDIR}/${t}_Ago3AS_SiwiS.pp8.q2.UA_VA.log"

	A=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.${t}.S.norm.bed.gz
	B=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.${t}.AS.norm.bed.gz
	jobname=${t}_Ago3S_SiwiAS.pp8.q2
	OUT=${OUTDIR}/${t}_Ago3S_SiwiAS
	[ ! -d ${OUT} ] && mkdir -p ${OUT}
	/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "${script} -i ${A} -j ${B} -n 2 -s bombyx -o ${OUT} -d ${indexFlag} -f normbed >${OUTDIR}/${t}_Ago3S_SiwiAS.pp8.q2.UA_VA.log"	
done
fi
[ $? == 0 ] && \
touch .status.${STEP}.pp8_UA_VA_cistrans
STEP=$((STEP+1))

declare -a PPPAIR=("Ago3AS_SiwiS" "Ago3S_SiwiAS")
declare -a TARGETS=("KNOWNTE")
OUTDIR=/home/wangw1/isilon_temp/BmN4/pp8_UA_VA_cistrans_summary
[ ! -d ${OUTDIR} ] && mkdir $OUTDIR

if [ ! -f .status.${STEP}.pp8_UA_VA_cistrans_summary ] 
then
for t in ${TARGETS[@]}
do
	for pp in ${PPPAIR[@]}
	do
		fn=/home/wangw1/isilon_temp/BmN4/pp8_UA_VA_cistrans/${t}_${pp}
		file=/home/wangw1/isilon_temp/BmN4/pp8_UA_VA_cistrans_summary/${t}_${pp}
		[ -f ${file}.pair.count.txt ] &&  rm ${file}.pair.count.txt
		for i in `ls ${fn}/*.UA_VA.zscore.out`
		do
		[ -f ${i} ] && /home/wangw1/git/smallRNA_analysis/g1Ut10A/UA_VA_rawscore.pl ${i} $OUTDIR >> ${file}.pair.count.txt
		done
		#sort -k1,1 -k2,2 -k3,3 -k4,4 $file.pair.count.txt | uniq >${file}.pair.count.uniq.txt
		#RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va ${file}.pair.count.uniq.txt ${t}_${pp} ${OUTDIR}
		[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va_color ${file}.pair.count.txt ${t}_${pp} ${OUTDIR}
	done
	
	
done
fi
[ $? == 0 ] && \
touch .status.${STEP}.pp8_UA_VA_cistrans_summary
STEP=$((STEP+1))
