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
fa=/home/xuj1/pipeline/common/fasta/dmel-all-chromosome-r5.5_TAS.fasta

STEP=1

ppscript=${PIPELINE_DIRECTORY}/g1Ut10A/pp8_SRA_vs_DEG_q2_ww1_zscore_sep_02272014_baseFraction.pl
sharedscript=${PIPELINE_DIRECTORY}/SRA_DEG_shared5end.pl
#use xkxh.norm.bed instead of transposon mappers
SRAINDIR=/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results/
DEGINDIR=/home/wangw1/isilon_temp/degradome/pipeline_output_12262013/

feature=$1
#FLY_PIRNA_CLUSTER
OUTDIR=/home/wangw1/isilon_temp/ipsmRNA/SRADEG_link/${feature}

[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}

echo -e "`date` "+$ISO_8601"\trun IP smallRNA and IP degradome association analysis..." >> $LOG
#declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "ago3Hets" "aubHets" "qinHets" "nosAgo3CDrescue" "nosAgo3WTrescue")
declare -a GT=("AubWTrescue" "AubCDrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "w1")
declare -a OX=("ox" "unox")
declare -a UNIQ=("uniq")

indexFlag=1

if [ ! -f ${OUTDIR}/.status.${STEP}.IPSRA_IPDEG ]
then 
	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
			for s in ${UNIQ[@]}
			do
																				 
				demapper2=${DEGINDIR}/Phil.DEG.AubIP.${t}.ovary.PE/bedIntersectWW/Phil.DEG.AubIP.${t}.ovary.PE.x_rRNA.dm3.sorted.f0x40.noS.5p.all.bed.ntm.collapse.${feature}.nta.mapper2.gz
				denormbed=${DEGINDIR}/Phil.DEG.AubIP.${t}.ovary.PE/bedIntersectWW/Phil.DEG.AubIP.${t}.${feature}.ovary.norm.bed.gz
								
				#change to norm.bed format for degradome
				#[ ! -s ${demapper2%.gz}.norm.bed.gz ] && ${PIPELINE_DIRECTORY}/mapper2gznormbed.pl ${demapper2} ${DEGINDIR}/Phil.DEG.AubIP.${t}.ovary.PE/bedIntersectWW/ && gzip ${demapper2%.gz}.norm.bed
				#[ ! -s ${denormbed} ] && mv ${demapper2%.gz}.norm.bed.gz ${denormbed}
							
				#Ago3IP piRNAs shared 5'end with AubIP degradome
				smnormbed=${SRAINDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz
#					
#				####script to check the shared 5' end between 5 script
#				[ -f ${smnormbed} ] && [ -f ${denormbed} ] && \
#				jobname=${t}_Ago3IPSRA${s}_${o}_AubIPDEG.5shared && \
#				jOUT=${OUTDIR}/${t}_Ago3IPSRA${s}_${o}_AubIPDEG_shared && \
#				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
#				/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "$sharedscript -i $smnormbed -j $denormbed -o ${jOUT} -f normbed"		
#				
#				###Ago3IP piRNAs complementary to AubIP degradome
#
#				####script to check the percent of piRNAs or degradome that are complementary to each other
				[ -f ${smnormbed} ] && [ -f ${denormbed} ] && \
				jobname=${t}_Ago3IPSRA${s}_${o}_AubIPDEG.PP && \
				jOUT=${OUTDIR}/${t}_Ago3IPSRA${s}_${o}_AubIPDEG_PP && \
				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
				/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "$ppscript -i ${smnormbed} -j ${denormbed} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed -a ${fa}"	

				
				#AubIP piRNAs shared 5'end with AubIP degradome
				smnormbed=${SRAINDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.norm.bed.gz
#				[ -f ${smnormbed} ] && [ -f ${denormbed} ] && \
#				jobname=${t}_AubIPSRA${s}_${o}_AubIPDEG.5shared && \
#				jOUT=${OUTDIR}/${t}_AubIPSRA${s}_${o}_AubIPDEG_shared && \
#				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
#				/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "$sharedscript -i $smnormbed -j $denormbed -o ${jOUT} -f normbed"	
				

				#AubIP piRNAs complementary to AubIP degradome
				[ -f ${smnormbed} ] && [ -f ${denormbed} ] && \
				jobname=${t}_AubIPSRA${s}_${o}_AubIPDEG.PP && \
				jOUT=${OUTDIR}/${t}_AubIPSRA${s}_${o}_AubIPDEG_PP && \
				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
				/home/wangw1/bin/submitsge 24 ${jobname} $OUTDIR "$ppscript -i ${smnormbed} -j ${denormbed} -n 2 -s fly -o ${jOUT} -d ${indexFlag} -f normbed -a ${fa}"	


			done
		done
	done
fi
[ $? == 0 ] && \
touch ${OUTDIR}/.status.${STEP}.IPSRA_IPDEG
STEP=$((STEP+1))







