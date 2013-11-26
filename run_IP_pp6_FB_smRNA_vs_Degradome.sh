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

script=${PIPELINE_DIRECTORY}/pp6_T_ww_len.pl


INDIR=$1 #this is the folder store all pipeline results outmost folders 
#/home/wangw1/isilon_temp/ipsmRNA/jia_pipeline_results
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log

STEP=1

declare -a GROUPGT=("AubIP_ago3cdwt_ox" "AubIP_ago3cdwt_unox" \
"AubIP_ago3cdw1_ox" "AubIP_ago3cdw1_unox" \
)

echo -e "`date` "+$ISO_8601"\trun pp8 UA and VA..." >> $LOG
declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "ago3Hets" "aubHets" "qinHets" "nosAgo3CDrescue" "nosAgo3WTrescue")
declare -a OX=("ox" "unox")
declare -a UNIQ=("uniq" "shared")

#uniq
#shared
indexFlag=1 #to indicate we need to build the index or not
OUTDIR=${OUT}/pp6_TOTAL
[ ! -f ${OUTDIR} ] && mkdir -p ${OUTDIR}
#touch ${OUT}/.status.${STEP}.pp6.SRA_vs_SRA
if [ ! -f ${OUT}/.status.${STEP}.pp6.SRA_vs_SRA ] 
then
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do
		for s in ${UNIQ[@]}
		do
			A=${INDIR}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.gz
			B=${INDIR}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}${s}_${o}_Ago3_Aub.pp6 && \
			jOUT=${OUTDIR}/${t}${s}_${o}_Ago3_Aub && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "${script} ${A} ${B} 2 ${jOUT} >${jOUT}/${t}${s}_${o}_Ago3_Aub.total.pp6.out"
		done
	done
done
#total
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do

			A=${INDIR}/Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts//Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.gz
			B=${INDIR}/Phil.SRA.AubIP.${t}.${o}.ovary.inserts/Phil.SRA.AubIP.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.gz
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}_${o}_Ago3_AubA.pp6 && \
			jOUT=${OUTDIR}/${t}_${o}_Ago3_Aub && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "${script} ${A} ${B} 2 ${jOUT} >${jOUT}/${t}_${o}_Ago3_Aub.total.pp6.out"
	done
done

fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.pp6.SRA_vs_SRA
STEP=$((STEP+1))

OUTDIR=${OUT}/FB
[ ! -f ${OUTDIR} ] && mkdir -p ${OUTDIR}
if [ ! -f ${OUT}/.status.${STEP}.FB ] 
then
	for i in `ls ${INDIR}/*.inserts/*.inserts.xkxh.transposon.mapper2.gz`
	do
		gt=${i##*/}
		gt=${gt/Phil.SRA./}
		gt=${gt/.ovary.*/}
		GTDIR=${OUTDIR}/${gt}
		[ ! -f ${GTDIR} ] && mkdir -p ${GTDIR}
		smmapper2=$i
		[ ! -f ${smmapper2%.gz}.23-29.roo.gz ] && \
		/home/wangw1/bin/submitsge 4 ${gt} $OUTDIR "${PIPELINE_DIRECTORY}/gzlenrangeselector.pl ${smmapper2} 23 29 >${smmapper2%.gz}.23-29 && \
			gzip ${smmapper2%.gz}.23-29 && \
			${PIPELINE_DIRECTORY}/FB.pl ${smmapper2%.gz}.23-29.gz ${GTDIR}"
	done
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.pp6_FB
STEP=$((STEP+1))

OUTDIR=${OUT}/pp6_FB
[ ! -f ${OUTDIR} ] && mkdir -p ${OUTDIR}
script=${PIPELINE_DIRECTORY}/pp6_T_ww_smRNA_vs_DEG.pl
TRN=("dodeca_satellite" "FBgn0000004_17.6" "FBgn0000005_297" "FBgn0000006_412" "FBgn0000007_1731" "FBgn0000155_roo" "FBgn0000199_blood" "FBgn0000224_BS" "FBgn0000349_copia" "FBgn0000481_Doc" "FBgn0000638_FB" "FBgn0000652_F_element" "FBgn0001100_G_element" "FBgn0001167_gypsy" "FBgn0001181_HB" "FBgn0001207_HMS_Beagle" "FBgn0001210_hobo" "FBgn0001249_I_element" "FBgn0001283_jockey" "FBgn0002697_mdg1" "FBgn0002698_mdg3" "FBgn0002745_micropia" "FBgn0002949_NOF" "FBgn0003007_opus" "FBgn0003122_pogo" "FBgn0003490_springer" "FBgn0003519_Stalker" "FBgn0003908_R1A1_element" "FBgn0003909_R2_element" "FBgn0004082_Tirant" "FBgn0004141_HeT_A" "FBgn0004904_TART_C_TART_B_TART_A" "FBgn0004905_S_element" "FBgn0005384_3S18" "FBgn0005673_1360" "FBgn0005773_Bari1" "FBgn0010103_aurora_element" "FBgn0010302_Burdock" "FBgn0014947_flea" "FBgn0014967_hopper" "FBgn0015786_Porto1" "FBgn0015945_GATE" "FBgn0020425_Helena" "FBgn0022937_Circe" "FBgn0023131_ZAM" "FBgn0026065_Idefix" "FBgn0026410_Tc1" "FBgn0026416_INE_1" "FBgn0040267_Transpac" "FBgn0041728_Rt1a" "FBgn0042231_X_element" "FBgn0042682_Rt1b" "FBgn0043055_Ivk" "FBgn0043969_diver" "FBgn0044355_Quasimodo" "FBgn0045970_Tabor" "FBgn0046110_Juan" "FBgn0061191_Tc3" "FBgn0061485_rover" "FBgn0061513_frogger" "FBgn0062343_Dm88" "FBgn0063369_transib4" "FBgn0063370_transib3" "FBgn0063371_transib2" "FBgn0063372_transib1" "FBgn0063394_rooA" "FBgn0063401_mariner2" "FBgn0063402_looper1" "FBgn0063425_jockey2" "FBgn0063426_invader5" "FBgn0063427_invader4" "FBgn0063428_invader3" "FBgn0063429_invader2" "FBgn0063430_invader1" "FBgn0063431_gypsy6" "FBgn0063432_gypsy5" "FBgn0063433_gypsy4" "FBgn0063434_gypsy3" "FBgn0063435_gypsy2" "FBgn0063436_gtwin" "FBgn0063439_diver2" "FBgn0063440_baggins" "FBgn0063447_accord" "FBgn0063450_Tom1" "FBgn0063454_Stalker3" "FBgn0063455_Stalker2" "FBgn0063466_S2" "FBgn0063467_Rt1c" "FBgn0063503_G6" "FBgn0063504_G5" "FBgn0063505_G4" "FBgn0063506_G3" "FBgn0063507_G2" "FBgn0063533_Doc3_element" "FBgn0063534_Doc2_element" "FBgn0063594_Cr1a" "FBgn0063755_Osvaldo" "FBgn0063782_accord2" "FBgn0063897_Stalker4" "FBgn0063900_Q_element" "FBgn0063917_McClintock" "FBgn0063919_Max_element" "FBgn0064134_Bari2" "FBgn0067381_hopper2" "FBgn0067382_gypsy9" "FBgn0067383_gypsy8" "FBgn0067384_gypsy7" "FBgn0067385_gypsy12" "FBgn0067386_gypsy11" "FBgn0067387_gypsy10" "FBgn0067405_R1_2" "FBgn0067418_Helitron" "FBgn0067419_G7" "FBgn0067420_Fw3" "FBgn0067421_Fw2" "FBgn0067624_BS3" "FBgn0069340_Tc1_2" "FBgn0069343_TAHRE" "FBgn0069433_G5A" "FBgn0069587_Doc4_element" "FBgnnnnnnnn_HMS_Beagle2" "mst40" "P-element" "stellate" "stellateHet" "suffix" "suste" "TARTA")
if [ ! -f ${OUT}/.status.${STEP}.pp6_FB_SRA_vs_SRA ] 
then
	for t in ${GT[@]}
	do
		for o in ${OX[@]}
		do
			for s in ${UNIQ[@]}				
			do
				jOUT=${OUTDIR}/${t}${s}_${o}_Ago3_Aub && \
				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
				paraFile=${jOUT}/${t}${s}_${o}_Ago3_Aub.parafile
				for fb in ${TRN[@]}
				do	
				
				A=${INDIR}/transposon_piRNA/FB/Ago3IP${s}.${t}.${o}/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.23-29.${fb}
				B=${INDIR}/transposon_piRNA/FB/AubIP${s}.${t}.${o}/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.23-29.${fb}
				[ -f ${A} ] && [ -f ${B} ] && \
				jobname=${t}${s}_${o}_Ago3_Aub.pp6 && \
				
				echo -ne "${script} ${A} ${B} 2 ${jOUT} ${fb} >>${jOUT}/${t}${s}_${o}_Ago3_Aub.FB.${fb}.pp6.out && ">>${paraFile} ####
				echo -e "${script} ${A} ${B} 2 ${jOUT} ${fb} >>${jOUT}/${t}${s}_${o}_Ago3_Aub.FB.pp6.out.temp">>${paraFile}
				
			done
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "ParaFly -c ${paraFile} -CPU 8 -failed_cmds ${paraFile}.failed_commands && \
			awk '{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${jOUT}/${t}${s}_${o}_Ago3_Aub.FB.pp6.out.temp >${jOUT}/${t}${s}_${o}_Ago3_Aub.FB.pp6.out"
		done
	done
done
#total
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do
				jOUT=${OUTDIR}/${t}_${o}_Ago3_Aub && \
				[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
				paraFile=${jOUT}/${t}_${o}_Ago3_Aub.parafile
				for fb in ${TRN[@]}
				do	
				
				A=${INDIR}/transposon_piRNA/FB/Ago3IP.${t}.${o}/Phil.SRA.Ago3IP.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.23-29.${fb}
				B=${INDIR}/transposon_piRNA/FB/AubIP.${t}.${o}/Phil.SRA.AubIP.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.23-29.${fb}
				[ -f ${A} ] && [ -f ${B} ] && \
				jobname=${t}_${o}_Ago3_Aub.pp6 && \
				
				echo -ne "${script} ${A} ${B} 2 ${jOUT} ${fb} >>${jOUT}/${t}_${o}_Ago3_Aub.FB.${fb}.pp6.out && ">>${paraFile}
				echo -e "${script} ${A} ${B} 2 ${jOUT} ${fb} >>${jOUT}/${t}_${o}_Ago3_Aub.FB.pp6.out.temp">>${paraFile}
				done
				/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "ParaFly -c ${paraFile} -CPU 8 -failed_cmds ${paraFile}.failed_commands && \
				awk '{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${jOUT}/${t}_${o}_Ago3_Aub.FB.pp6.out.temp > ${jOUT}/${t}_${o}_Ago3_Aub.FB.pp6.out"
				

	done
done
				
fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.pp6_FB_SRA_vs_SRA
STEP=$((STEP+1))	


#use Chi-square enriched species 
declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue")
declare -a OX=("ox" "unox")
declare -a UNIQ=("enriched")
OUTDIR=${OUT}/pp6_TOTAL_enriched
[ ! -f ${OUTDIR} ] && mkdir -p ${OUTDIR}
#touch ${OUT}/.status.${STEP}.pp6.SRA_vs_SRA
if [ ! -f ${OUT}/.status.${STEP}.pp6.enriched.SRA_vs_SRA ] 
then
for t in ${GT[@]}
do
	for o in ${OX[@]}
	do
		for s in ${UNIQ[@]}
		do
			A=/home/wangw1/isilon_temp/ipsmRNA/chi-square_enriched/Phil.SRA.Ago3IP${s}.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.23-29
			B=/home/wangw1/isilon_temp/ipsmRNA/chi-square_enriched/Phil.SRA.AubIP${s}.${t}.${o}.ovary.inserts.xkxh.transposon.mapper2.23-29
			[ -f ${A} ] && [ -f ${B} ] && \
			jobname=${t}${s}_${o}_Ago3_Aub.pp6 && \
			jOUT=${OUTDIR}/${t}${s}_${o}_Ago3_Aub && \
			[ ! -d ${jOUT} ] && mkdir -p ${jOUT} && \
			/home/wangw1/bin/submitsge 8 ${jobname} $OUTDIR "${script} ${A} ${B} 2 ${jOUT} >${jOUT}/${t}${s}_${o}_Ago3_Aub.total.pp6.out"
		done
	done
done


fi
[ $? == 0 ] && \
touch ${OUT}/.status.${STEP}.pp6.enriched.SRA_vs_SRA
STEP=$((STEP+1))

