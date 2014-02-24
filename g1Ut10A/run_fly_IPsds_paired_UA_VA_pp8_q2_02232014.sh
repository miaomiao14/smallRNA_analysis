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

indexFlag=$1 #this is the folder store all pipeline results outmost folders 
#/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results


STEP=1

script=${PIPELINE_DIRECTORY}/g1Ut10A/pp8_q2_ww1_zscore_sep_02232014.pl

#run UA_VA analysis
#separate sense and antisense from normbed
#Phil.AubIP.AubCDrescue.unox.ovary.inserts.xkxh.norm.bed.transposons.gz

#uniq bound piRNAs exlude Piwi bound(less piRNA input, to be consistent with the analysis done in 2012)
#use xkxh.norm.bed instead of transposon mappers
INDIR=/home/wangw1/nearline/mpo/IP
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log
declare -a GT=("ago3Hets" "aubHets" "qinHets" "wt")
declare -a OX=("unox")
declare -a UNIQ=("uniq")
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results/transposon_piRNA/UA_VA_cistrans_allpiRNAs_exclude_Piwi
[ ! -d ${OUTDIR} ] && mkdir $OUTDIR
#touch ${OUT}/.status.${STEP}.all_piRNA_exclude_Piwi.UA_VA
if [ ! -f ${OUT}/.status.${STEP}.all_piRNA_exclude_Piwi.UA_VA_cistrans ]
then 
	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
			for s in ${UNIQ[@]}
			do
				A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				[ -f ${A} ] && [ -f ${B} ] && \
				jobname=${t}${s}_${o}_Ago3_Aub.pp8.q2 && \
				jOUT=${OUTDIR}/${t}${s}_${o}_Ago3_Aub && \
				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
				/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed >${jOUT}/${t}${s}_${o}_Ago3_Aub.pp8.q2.UA_VA.log"
#			
				done
		done
	done

fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.all_piRNA_exclude_Piwi.UA_VA_cistrans
STEP=$((STEP+1))


declare -a PPPAIR=("Ago3_Aub")
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results/transposon_piRNA/UA_VA_allpiRNAs_exclude_Piwi
SUMMARYOUTDIR=/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results/transposon_piRNA/UA_VA_allpiRNAs_exclude_Piwi_summary
[ ! -d ${SUMMARYOUTDIR} ] && mkdir $SUMMARYOUTDIR
#touch .status.${STEP}.pp8_all_piRNA_exclude_Piwi.UA_VA_summary
if [ ! -f .status.${STEP}.pp8_all_piRNA_exclude_Piwi.UA_VA_summary ] 
then
	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
			for s in ${UNIQ[@]}
			do
				for pp in ${PPPAIR[@]}
				do
					fn=${OUTDIR}/${t}${s}_${o}_${pp}
					file=${SUMMARYOUTDIR}/${t}${s}_${o}_${pp}
					[ -f ${file}.pair.count.txt ] &&  rm ${file}.pair.count.txt 
					for i in `ls ${fn}/*.VA.pp`
					do
						filename=${i##*/}
					pairname=`basename ${filename} .VA.pp`
					arrName=(${pairname//\./ })
					pairnameNew=${arrName[1]}"-"${arrName[0]}
				
					[ -f $i ] && awk -v gt=${pairnameNew} '{OFS="\t"}{print gt,$0}' ${i} >${i}.gt && \
					/home/wangw1/git/smallRNA_analysis/Ping_Pong/UA_VA_rawscore_from_ppscore.pl ${i}.gt $OUTDIR >> ${file}.pair.count.txt
					#sort -k1,1 -k2,2 -k3,3 -k4,4 $file.pair.count.txt | uniq >${file}.pair.count.uniq.txt
					#[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va ${file}.pair.count.txt ${t}${s}_${o}_${pp} ${SUMMARYOUTDIR}
					[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va_from_ppscore_color ${file}.pair.count.txt ${t}${s}_${o}_${pp} ${SUMMARYOUTDIR}
					done
				done
			done
		done		
	done

fi
[ $? == 0 ] && \
touch .status.${STEP}.pp8_all_piRNA_exclude_Piwi.UA_VA_summary
STEP=$((STEP+1))



#IPed piRNA abundance normalzation



