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
declare -a GT=("ago3Hets" "aubHets" "qinHets")
declare -a OX=("unox")
declare -a UNIQ=("uniq")
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cistrans_allpiRNAs_exclude_Piwi
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


#run UA_VA analysis
#separate sense and antisense from normbed
#Phil.AubIP.AubCDrescue.unox.ovary.inserts.xkxh.norm.bed.transposons.gz
echo -e "`date` "+$ISO_8601"\t1U10A and 1V10A analysis..." >> $LOG
echo -e "`date` "+$ISO_8601"\tseparate the sense and antisense transposon mappers in norm.bed format..." >> $LOG
#touch ${OUT}/.status.transposon_piRNA.S_AS
INDIR=/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans
[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}
[ ! -f ${OUT}/.status.transposon_piRNA.S_AS ] && \
paraFile=${OUTDIR}/UA_VA_cis_trans.${RANDOM}.para && \
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
STEP=$((STEP+1))

echo -e "`date` "+$ISO_8601"\trun pp8 UA and VA..." >> $LOG
declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "ago3Hets" "aubHets" "qinHets" "nosAgo3CDrescue" "nosAgo3WTrescue")
declare -a OX=("ox" "unox")
declare -a UNIQ=("uniq")

#uniq
#shared
indexFlag=1 #to indicate we need to build the index or not
	
#touch ${OUT}/.status.${STEP}.transposon_piRNA.UA_VA_cis_trans
if [ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.UA_VA_cis_trans ] 
then
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do
		for s in ${UNIQ[@]}
		do
			A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.sense.gz
			B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.antisense.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}${s}_${o}_Ago3S_AubAS.pp8.q2 && \
			jOUT=${OUTDIR}/${t}${s}_${o}_Ago3S_AubAS && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed >${jOUT}/${t}${s}_${o}_Ago3S_AubAS.pp8.q2.UA_VA_cis_trans.log"
#			
			A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.antisense.gz
			B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.sense.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}${s}_${o}_Ago3AS_AubS.pp8.q2 && \
			jOUT=${OUTDIR}/${t}${s}_${o}_Ago3AS_AubS && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed >${jOUT}/${t}${s}_${o}_Ago3AS_AubS.pp8.q2.UA_VA_cis_trans.log"			
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
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed >${jOUT}/${t}_${o}_Ago3S_AubAS.pp8.q2.UA_VA_cis_trans.log"
			
			A=${INDIR}/Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.antisense.gz
			B=${INDIR}/Phil.SRA.AubIP.${t}.${o}.ovary.inserts/Phil.SRA.AubIP.${t}.${o}.ovary.inserts.xkxh.norm.bed.transposons.sense.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}_${o}_Ago3AS_AubS.pp8.q2 && \
			jOUT=${OUTDIR}/${t}_${o}_Ago3AS_AubS && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed >${jOUT}/${t}_${o}_Ago3AS_AubS.pp8.q2.UA_VA_cis_trans.log"			

	done
done
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.transposon_piRNA.UA_VA_cis_trans
STEP=$((STEP+1))


declare -a PPPAIR=("Ago3AS_AubS" "Ago3S_AubAS")

OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans
SUMMARYOUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans_summary
[ ! -d ${SUMMARYOUTDIR} ] && mkdir $SUMMARYOUTDIR
touch ${OUT}/.status.${STEP}.pp8_UA_VA_cis_trans_summary
if [ ! -f ${OUT}/.status.${STEP}.pp8_UA_VA_cis_trans_summary ] 
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
				for i in `ls ${fn}/*.UA_VA.zscore.out`
				do
				#filename=${i##*/}
				#pairname=`basename ${filename} .VA.pp`
				[ -f $i ] && /home/wangw1/git/smallRNA_analysis/Ping_Pong/UA_VA_rawscore.pl ${i} $OUTDIR >> ${file}.pair.count.txt 
				done
				#sort -k1,1 -k2,2 -k3,3 -k4,4 $file.pair.count.txt | uniq >${file}.pair.count.uniq.txt
				#[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va ${file}.pair.count.txt ${t}${s}_${o}_${pp} ${SUMMARYOUTDIR}
				[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va_color ${file}.pair.count.txt ${t}${s}_${o}_${pp} ${SUMMARYOUTDIR}
			done
		done
	done		
done

for t in ${GT[@]}
do
	for o in ${OX[@]}
	do

			for pp in ${PPPAIR[@]}
			do
				fn=${OUTDIR}/${t}_${o}_${pp} 
				file=${SUMMARYOUTDIR}/${t}_${o}_${pp} 
				[ -f ${file}.pair.count.txt ] &&  rm ${file}.pair.count.txt 
				for i in `ls ${fn}/*.UA_VA.zscore.out`
				do
					[ -f $i ] && /home/wangw1/git/smallRNA_analysis/Ping_Pong/UA_VA_rawscore.pl ${i} $OUTDIR >> ${file}.pair.count.txt 
				done
				#sort -k1,1 -k2,2 -k3,3 -k4,4 $file.pair.count.txt | uniq >${file}.pair.count.uniq.txt
				#[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va ${file}.pair.count.txt ${t}_${o}_${pp} ${SUMMARYOUTDIR}
				[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va_color ${file}.pair.count.txt ${t}_${o}_${pp} ${SUMMARYOUTDIR} 
			done

	done		
done
fi


[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.pp8_UA_VA_cis_trans_summary
STEP=$((STEP+1))


#use xkxh.norm.bed instead of transposon mappers
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans_allpiRNAs
#touch ${OUT}/.status.${STEP}.all_piRNA.UA_VA_cis_trans
if [ ! -f ${OUT}/.status.${STEP}.all_piRNA.UA_VA_cis_trans ]
then 
	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
			for s in ${UNIQ[@]}
			do
				A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				[ -f ${A} ] && [ -f ${B} ] && \
				jobname=${t}${s}_${o}_Ago3_Aub.pp8.q2 && \
				jOUT=${OUTDIR}/${t}${s}_${o}_Ago3_Aub && \
				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
				/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed >${jOUT}/${t}${s}_${o}_Ago3_Aub.pp8.q2.UA_VA_cis_trans.log"
#			
				done
		done
	done
#total

	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
	
				A=${INDIR}/Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				B=${INDIR}/Phil.SRA.AubIP.${t}.${o}.ovary.inserts/Phil.SRA.AubIP.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz 
				[ -f ${A} ] && [ -f ${B} ] && \
				jobname=${t}_${o}_Ago3_Aub.pp8.q2 && \
				jOUT=${OUTDIR}/${t}_${o}_Ago3_Aub && \
				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
				/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "$script -i ${A} -j ${B} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed >${jOUT}/${t}_${o}_Ago3_Aub.pp8.q2.UA_VA_cis_trans.log"
	
		done
	done
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.all_piRNA.UA_VA_cis_trans
STEP=$((STEP+1))


declare -a PPPAIR=("Ago3_Aub")
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans_allpiRNAs
SUMMARYOUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans_allpiRNAs_summary
[ ! -d ${SUMMARYOUTDIR} ] && mkdir $SUMMARYOUTDIR
if [ ! -f .status.${STEP}.pp8_all_piRNA.UA_VA_cis_trans_summary ] 
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
					for i in `ls ${fn}/*.UA_VA.zscore.out`
					do
					[ -f $i ] && /home/wangw1/git/smallRNA_analysis/Ping_Pong/UA_VA_rawscore.pl ${i} $OUTDIR >> ${file}.pair.count.txt 
					done
					#sort -k1,1 -k2,2 -k3,3 -k4,4 $file.pair.count.txt | uniq >${file}.pair.count.uniq.txt
					#[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va ${file}.pair.count.txt ${t}${s}_${o}_${pp} ${SUMMARYOUTDIR}
					[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va_color ${file}.pair.count.txt ${t}${s}_${o}_${pp} ${SUMMARYOUTDIR}
	
				done
			done
		done		
	done
	
	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
	
				for pp in ${PPPAIR[@]}
				do
					fn=${OUTDIR}/${t}_${o}_${pp}
					file=${SUMMARYOUTDIR}/${t}_${o}_${pp}
					[ -f ${file}.pair.count.txt ] &&  rm ${file}.pair.count.txt 
					for i in `ls ${fn}/*.UA_VA.zscore.out`
					do
						[ -f $i ] && /home/wangw1/git/smallRNA_analysis/Ping_Pong/UA_VA_rawscore.pl ${i} $OUTDIR >> ${file}.pair.count.txt 
					done
					#sort -k1,1 -k2,2 -k3,3 -k4,4 $file.pair.count.txt | uniq >${file}.pair.count.uniq.txt
					#[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va ${file}.pair.count.txt ${t}_${o}_${pp} ${SUMMARYOUTDIR}
					[ -f ${file}.pair.count.txt ] && RRR /home/wangw1/git/smallRNA_analysis/Ping_Pong/ping_pong_plots.r plot_ua_va_color ${file}.pair.count.txt ${t}_${o}_${pp} ${SUMMARYOUTDIR} 
				done
	
		done		
	done
fi
[ $? == 0 ] && \
touch .status.${STEP}.pp8_all_piRNA.UA_VA_cis_trans_summary
STEP=$((STEP+1))


declare -a PPPAIR=("Ago3_Aub")
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans_allpiRNAs_exclude_Piwi
SUMMARYOUTDIR=/home/wangw1/isilon_temp/ipsmRNA/transposon_piRNA/UA_VA_cis_trans_allpiRNAs_exclude_Piwi_summary
[ ! -d ${SUMMARYOUTDIR} ] && mkdir $SUMMARYOUTDIR
#touch .status.${STEP}.pp8_all_piRNA_exclude_Piwi.UA_VA_cis_trans_summary
if [ ! -f .status.${STEP}.pp8_all_piRNA_exclude_Piwi.UA_VA_cis_trans_summary ] 
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
touch .status.${STEP}.pp8_all_piRNA_exclude_Piwi.UA_VA_cis_trans_summary
STEP=$((STEP+1))



#IPed piRNA abundance normalzation



