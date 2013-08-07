#!/bin/bash


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
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log

declare -a NORMFACTORTYPE=("nnc" "seqDep")

STEP=1
OUTDIR1=${INDIR}/transposon_piRNA/polar
[ ! -d $OUTDIR1 ] && mkdir -p ${OUTDIR1}
echo -e "`date` "+$ISO_8601"\tDraw polar histogram of sense fraction" >> $LOG
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.senseFraction ] && \
for i in `ls ${INDIR}/*.inserts/output/*.transposon.list`
do 

	FILE=${i##*/}
	insertsname=`basename $FILE .transposon.list`

#polarHistogram for sense fraction
${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/polarHistogram_sense_fraction.r plot_PolarHisto_senseFraction $i ${OUTDIR1}
done
[ $? == 0 ] && \	
	touch ${OUT}/.status.${STEP}.transposon_piRNA.senseFraction
STEP=$((STEP+1))


echo -e "`date` "+$ISO_8601"\tDraw length distribution of transposon piRNAs" >> $LOG
OUTDIR2=${INDIR}/transposon_piRNA/lendis
[ ! -d $OUTDIR2 ] && mkdir -p ${OUTDIR2}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.lendis2 ] && \
paraFile=${OUTDIR2}/${RANDOM}.drawlendis2.para && \
for i in `ls ${INDIR}/*.inserts/*.xkxh.transposon.mapper2.gz`
do 	
	FILE=${i##*/}
	insertsname=`basename $FILE .xkxh.transposon.mapper2.gz`
	inserts=${FILE%%.inserts.*}
	inserts=${inserts}.inserts
	nfnnc=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	nfdep=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	
	declare -a NORMFACTOR=(${nfnnc} ${nfdep}) #not in use now
	
	echo -ne "${PIPELINE_DIRECTORY}/lendis2.pl ${i} $OUTDIR2 nnc $nfnnc &&" >>${paraFile}
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_lendis2 ${OUTDIR2}/$insertsname.xkxh.transposon.mapper2.nnc.lendis2 ${insertsname}" >>${paraFile}	
	echo -ne "${PIPELINE_DIRECTORY}/lendis2.pl ${i} $OUTDIR2 seqDep $nfnnc &&" >>${paraFile}
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_lendis2 ${OUTDIR2}/$insertsname.xkxh.transposon.mapper2.seqDep.lendis2 ${insertsname}" >>${paraFile}		
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA.lendis2
STEP=$((STEP+1))



#declare -a gt_g1_ox=("aubvasAgo3CDrescue.ox" "ago3MutsWW.ox" "aubvasAgo3WTrescue.ox")
#declare -a gt_g1_unox=("aubvasAgo3CDrescue.unox" "ago3MutsWW.unox" "aubvasAgo3WTrescue.unox")
#declare -a gt_g2_ox=("AubCDrescue.ox" "AubMutsWW.ox" "AubWTrescue.ox")
#declare -a gt_g2_unox=("AubCDrescue.unox" "AubMutsWW.unox" "AubWTrescue.unox")
#
#declare -a gt_cor1_ox=("ago3MutsWW.ox" "ago3MutsCJ.ox")
#declare -a gt_cor1_unox=("ago3MutsWW.unox" "ago3MutsCJ.unox")
#
#declare -a gt_cor2_ox=("nosAgo3CDrescue.ox" "aubvasAgo3CDrescue.ox")
#declare -a gt_cor2_unox=("nosAgo3CDrescue.unox" "aubvasAgo3CDrescue.unox")

declare -a GROUPGT=("pago3_ox" "pago3_unox" "ago3_ox" "ago3_unox" "aub_ox" "aub_unox" "ago3mut_cor1_ox" "ago3mut_cor1_unox" "ago3CD_cor2_ox" "ago3CD_cor2_unox")
declare -a pago3_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.nosAgo3WTrescue.ox.ovary.inserts")
declare -a pago3_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.nosAgo3WTrescue.unox.ovary.inserts")
declare -a ago3_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts")
declare -a ago3_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts")
declare -a aub_ox=("Phil.SRA.AubCDrescue.ox.ovary.inserts" "Phil.SRA.AubMutsWW.ox.ovary.inserts" "Phil.SRA.AubWTrescue.ox.ovary.inserts")
declare -a aub_unox=("Phil.SRA.AubCDrescue.unox.ovary.inserts" "Phil.SRA.AubMutsWW.unox.ovary.inserts" "Phil.SRA.AubWTrescue.unox.ovary.inserts")



declare -a ago3mut_cor1_ox=("Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.ago3MutsCJ.ox.ovary.inserts")
declare -a ago3mut_cor1_unox=("Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.ago3MutsCJ.unox.ovary.inserts")

declare -a ago3CD_cor2_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts")
declare -a ago3CD_cor2_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts")

echo -e "`date` "+$ISO_8601"\tDraw paired length distribution of transposon piRNAs" >> $LOG
OUTDIR3=${INDIR}/transposon_piRNA/paired_lendis
[ ! -d $OUTDIR3 ] && mkdir -p ${OUTDIR3}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.paired.lendis2 ] && \
paraFile=${OUTDIR3}/${RANDOM}.drawpairedlendis2.para && \
for g in "${GROUPGT[@]}"
do
	SUBGROUP="$g[@]"
	echo $g  >> $LOG
	#echo ${#(!SUBGROUP)} >> $LOG #not the value
	#declare -a MAPPER2NNCLENDIS=()
	#declare -a MAPPER2UNIQLENDIS=()
	#MAPPER2NNCLENDIS=${MAPPER2NNCLENDIS}","${OUTDIR2}/${t}.xkxh.transposon.mapper2.nnc.lendis2 
	#MAPPER2UNIQLENDIS=${MAPPER2UNIQLENDIS}","${OUTDIR2}/${t}.uniqmap.xkxh.transposon.mapper2.nnc.lendis2 
	
	
	
	for NF in "${NORMFACTORTYPE[@]}"
	do
		lendisFile=${OUTDIR3}/${g}.${NF}.${RANDOM}.lendis2
		count=1
		for t in ${!SUBGROUP}
		do
		cat ${OUTDIR2}/${t}.xkxh.transposon.mapper2.${NF}.lendis2| awk -v gt=$t -v rank=$count '{OFS="\t"}{print gt,$1,$2,$3,rank}' >> ${lendisFile}
		count=$(($count+1))	 
		done
		echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_paired_lendis2 ${lendisFile} ${OUTDIR3}" >>${paraFile}
	done
	
	
	
	for NF in "${NORMFACTORTYPE[@]}"
	do
		lendisFile=${OUTDIR3}/${g}.${NF}.uniqmap.${RANDOM}.lendis2
		count=1
		for t in ${!SUBGROUP}
		do
		cat ${OUTDIR2}/${t}.uniqmap.xkxh.transposon.mapper2.${NF}.lendis2| awk -v gt=$t -v rank=$count '{OFS="\t"}{ print gt".uniqmap",$1,$2,$3,rank}' >> ${lendisFile}
		count=$(($count+1))		 
		done
		echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_paired_lendis2 ${lendisFile} ${OUTDIR3}" >>${paraFile}
	done
	
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA.paired.lendis2
STEP=$((STEP+1))



echo -e "`date` "+$ISO_8601"\tDraw transposon piRNA abundance and zscore barplot..." >> $LOG
OUTDIR4=${INDIR}/transposon_piRNA/abundance_zscore
[ ! -d $OUTDIR4 ] && mkdir -p ${OUTDIR4}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.abundance_zscore ] && \
paraFile=${OUTDIR4}/${RANDOM}.drawbarplotAbundanceZscore.para && \
for i in `ls ${INDIR}/*.inserts/*.xkxh.transposon.mapper2.gz`
do 	
	FILE=${i##*/}
	insertsname=`basename $FILE .xkxh.transposon.mapper2.gz`
	inserts=${FILE%%.inserts.*}
	inserts=${inserts}.inserts
	nnc=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	seqDep=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	
	declare -a NORMFACTOR=(${nfnnc} ${nfdep}) #not in use now
	
	totalZscore=`cat ${INDIR}/pp6_FB/${insertsname}/$insertsname.total.pp6.out|cut -f2`
	TZ=`printf "%0.2f" $totalZscore`
	#declare -a NORMFACTORTYPE=("nnc" "seqDep")	
	for NF in "${NORMFACTORTYPE[@]}"
	do
		echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_transposon_abundance_zscore_barplot ${INDIR}/${inserts}/output/${insertsname}.transposon.list ${INDIR}/pp6_FB/${insertsname}/$insertsname.FB.pp6.out $TZ ${!NF} $NF $OUTDIR4" >>${paraFile}	
	done
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA.abundance_zscore
STEP=$((STEP+1))




OUTDIR5=${INDIR}/transposon_piRNA/polar_allTrn
[ ! -d $OUTDIR5 ] && mkdir -p ${OUTDIR5}
echo -e "`date` "+$ISO_8601"\tDraw polar histogram of sense fraction(including transposon families not in groups)" >> $LOG
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.senseFractionallTrn ] && \
for i in `ls ${INDIR}/*.inserts/output/*.transposon.list`
do 

	FILE=${i##*/}
	insertsname=`basename $FILE .transposon.list`
${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/polarHistogram_sense_fraction.r plot_PolarHisto_senseFraction_allTrn $i ${OUTDIR5}
done
[ $? == 0 ] && \	
	touch ${OUT}/.status.${STEP}.transposon_piRNA.senseFractionallTrn
STEP=$((STEP+1))



declare -a GROUPGT=("pago3_ox" "pago3_unox" "ago3_ox" "ago3_unox" "aub_ox" "aub_unox" "ago3mut_cor1_ox" "ago3mut_cor1_unox" "ago3CD_cor2_ox" "ago3CD_cor2_unox")
declare -a pago3_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.nosAgo3WTrescue.ox.ovary.inserts")
declare -a pago3_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.nosAgo3WTrescue.unox.ovary.inserts")
declare -a ago3_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts")
declare -a ago3_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts")
declare -a aub_ox=("Phil.SRA.AubCDrescue.ox.ovary.inserts" "Phil.SRA.AubMutsWW.ox.ovary.inserts" "Phil.SRA.AubWTrescue.ox.ovary.inserts")
declare -a aub_unox=("Phil.SRA.AubCDrescue.unox.ovary.inserts" "Phil.SRA.AubMutsWW.unox.ovary.inserts" "Phil.SRA.AubWTrescue.unox.ovary.inserts")



declare -a ago3mut_cor1_ox=("Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.ago3MutsCJ.ox.ovary.inserts")
declare -a ago3mut_cor1_unox=("Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.ago3MutsCJ.unox.ovary.inserts")

declare -a ago3CD_cor2_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts")
declare -a ago3CD_cor2_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts")


INDIR6=${INDIR}/circos
OUTDIR6=/home/wangw1/src/circos-0.56/fly
OUTDIR6_2=${INDIR}/transposon_piRNA/circos_conf
#/home/wangw1/src/circos-0.56/bin/circos -conf /home/wangw1/src/circos-0.56/fly/etc/
BINSIZE=$2
declare -a NORMFACTORTYPE=("nnc" "seqDep")
declare -a METRICTYPE=("max" "mean")
[ ! -d $OUTDIR6/circosPlot ] && mkdir -p ${OUTDIR6}/circosPlot
[ ! -d ${OUTDIR6_2} ] && mkdir -p ${OUTDIR6_2}
echo -e "`date` "+$ISO_8601"\tCreate circos configuration files and generate circos plot" >> $LOG
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.circosplot ] && \
paraFile=${OUTDIR6_2}/${RANDOM}.drawcircosplot.para && \
for g in "${GROUPGT[@]}"
do
	SUBGROUP="$g[@]"
	
	for NF in ${NORMFACTORTYPE[@]}
	do 
		for MT in ${METRICTYPE[@]}
		do
			declare -a TARGETS=() 
			declare -a TARGETSUNIQ=() 
			for t in ${!SUBGROUP}
			do			
			TARGETS=${TARGETS}" "${INDIR6}/$t.${BINSIZE}/${t}.*.${NF}.antisense.${MT}.${BINSIZE}.circos.txt
			TARGETS=${TARGETS}" "${INDIR6}/$t.${BINSIZE}/${t}.*.${NF}.sense.${MT}.${BINSIZE}.circos.txt
			TARGETSUNIQ=${TARGETSUNIQ}" "${INDIR6}/$t.uniqmap.${BINSIZE}/${t}.uniqmap.*.${NF}.antisense.${MT}.${BINSIZE}.circos.txt
			TARGETSUNIQ=${TARGETSUNIQ}" "${INDIR6}/$t.uniqmap.${BINSIZE}/${t}.uniqmap.*.${NF}.sense.${MT}.${BINSIZE}.circos.txt
			done
			
			#echo ${TARGETS[@]} >> $LOG
			#echo ${TARGETSUNIQ[@]} >> $LOG
			#### can not use paraFly, because circos plot produce the plot with the same name: circos.svg and circos.png
			#### multithreading will cause conflict
			
#			echo -ne " ${PIPELINE_DIRECTORY}/circos_conf_gnr.sh ${TARGETS[@]} ${g}.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			#echo -ne " mv ${OUTDIR6}/etc/file.conf ${OUTDIR6}/etc/${g}.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			echo -ne " cd ${OUTDIR6} && " >>${paraFile}
#			echo -ne " /home/wangw1/src/circos-0.56/bin/circos -conf ./etc/${g}.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			echo -ne " mv ${OUTDIR6}/circos.svg ${OUTDIR6}/circosPlot/${g}.${NF}.${MT}.${BINSIZE}.svg && " >>${paraFile}
#			echo -e " mv ${OUTDIR6}/circos.png ${OUTDIR6}/circosPlot/${g}.${NF}.${MT}.${BINSIZE}.png " >>${paraFile}
#			
#			echo -ne " ${PIPELINE_DIRECTORY}/circos_conf_gnr.sh ${TARGETSUNIQ[@]} ${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			#echo -ne " mv ${OUTDIR6}/etc/file.conf ${OUTDIR6}/etc/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			echo -ne " cd ${OUTDIR6} && " >>${paraFile}
#			echo -ne " /home/wangw1/src/circos-0.56/bin/circos -conf ./etc/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			echo -ne " mv ${OUTDIR6}/circos.svg ${OUTDIR6}/circosPlot/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.svg && " >>${paraFile}
#			echo -e " mv ${OUTDIR6}/circos.png ${OUTDIR6}/circosPlot/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.png " >>${paraFile}
			
			${PIPELINE_DIRECTORY}/circos_conf_gnr.sh ${TARGETS[@]} ${g}.${NF}.${MT}.${BINSIZE}.conf && \

			cd ${OUTDIR6} && \
			/home/wangw1/src/circos-0.56/bin/circos -conf ./etc/${g}.${NF}.${MT}.${BINSIZE}.conf && \
			mv ${OUTDIR6}/circos.svg ${OUTDIR6}/circosPlot/${g}.${NF}.${MT}.${BINSIZE}.svg && \
			mv ${OUTDIR6}/circos.png ${OUTDIR6}/circosPlot/${g}.${NF}.${MT}.${BINSIZE}.png
			
			${PIPELINE_DIRECTORY}/circos_conf_gnr.sh ${TARGETSUNIQ[@]} ${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && \
			#echo -ne " mv ${OUTDIR6}/etc/file.conf ${OUTDIR6}/etc/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
			cd ${OUTDIR6} && \
			/home/wangw1/src/circos-0.56/bin/circos -conf ./etc/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && \
			mv ${OUTDIR6}/circos.svg ${OUTDIR6}/circosPlot/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.svg && \
			mv ${OUTDIR6}/circos.png ${OUTDIR6}/circosPlot/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.png
			
			
		done
	done
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA.circosplot
STEP=$((STEP+1))


echo -e "`date` "+$ISO_8601"\tDraw length distribution of transposon piRNAs, excluding roo transposon families" >> $LOG
OUTDIR7=${INDIR}/transposon_piRNA/lendisWOroo
[ ! -d $OUTDIR7 ] && mkdir -p ${OUTDIR7}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA_WOroo.lendis2 ] && \
paraFile=${OUTDIR7}/${RANDOM}.drawlendis2WOroo.para && \
for i in `ls ${INDIR}/*.inserts/*.xkxh.transposon.mapper2.gz`
do 	
	FILE=${i##*/}
	insertsname=`basename $FILE .xkxh.transposon.mapper2.gz`
	inserts=${FILE%%.inserts.*}
	inserts=${inserts}.inserts
	nfnnc=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	nfdep=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	
	declare -a NORMFACTOR=(${nfnnc} ${nfdep}) #not in use now
	
	echo -ne "[ ! -f ${OUTDIR7}/$insertsname.xkxh.transposon.mapper2.nnc.lendis2 ] && ${PIPELINE_DIRECTORY}/lendis2_WOroo.pl ${i} $OUTDIR7 nnc $nfnnc &&" >>${paraFile}
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_lendis2 ${OUTDIR7}/$insertsname.xkxh.transposon.mapper2.nnc.lendis2 ${insertsname}" >>${paraFile}	
	echo -ne "[ ! -f ${OUTDIR7}/$insertsname.xkxh.transposon.mapper2.seqDep.lendis2 ] && ${PIPELINE_DIRECTORY}/lendis2_WOroo.pl ${i} $OUTDIR7 seqDep $nfnnc &&" >>${paraFile}
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_lendis2 ${OUTDIR7}/$insertsname.xkxh.transposon.mapper2.seqDep.lendis2 ${insertsname}" >>${paraFile}		
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 4 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA_WOroo.lendis2
STEP=$((STEP+1))


echo -e "`date` "+$ISO_8601"\tDraw paired length distribution of transposon piRNAs WOroo" >> $LOG
OUTDIR8=${INDIR}/transposon_piRNA/paired_lendis_WOroo
[ ! -d $OUTDIR8 ] && mkdir -p ${OUTDIR8}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.paired.lendis2WOroo ] && \
paraFile=${OUTDIR8}/${RANDOM}.drawpairedlendis2WOroo.para && \
for g in "${GROUPGT[@]}"
do
	SUBGROUP="$g[@]"
	echo $g  >> $LOG
	#echo ${#(!SUBGROUP)} >> $LOG #not the value
	#declare -a MAPPER2NNCLENDIS=()
	#declare -a MAPPER2UNIQLENDIS=()
	#MAPPER2NNCLENDIS=${MAPPER2NNCLENDIS}","${OUTDIR2}/${t}.xkxh.transposon.mapper2.nnc.lendis2 
	#MAPPER2UNIQLENDIS=${MAPPER2UNIQLENDIS}","${OUTDIR2}/${t}.uniqmap.xkxh.transposon.mapper2.nnc.lendis2 
	
	
	
	for NF in "${NORMFACTORTYPE[@]}"
	do
		lendisFile=${OUTDIR8}/${g}.${NF}.${RANDOM}.WOroo.lendis2
		count=1
		for t in ${!SUBGROUP}
		do
		cat ${OUTDIR7}/${t}.xkxh.transposon.mapper2.${NF}.lendis2| awk -v gt=$t -v rank=$count '{OFS="\t"}{print gt,$1,$2,$3,rank}' >> ${lendisFile}
		count=$(($count+1))	 
		done
		echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_paired_lendis2 ${lendisFile} ${OUTDIR8}" >>${paraFile}
	done
	
	
	
	for NF in "${NORMFACTORTYPE[@]}"
	do
		lendisFile=${OUTDIR8}/${g}.${NF}.uniqmap.${RANDOM}.WOroo.lendis2
		count=1
		for t in ${!SUBGROUP}
		do
		cat ${OUTDIR7}/${t}.uniqmap.xkxh.transposon.mapper2.${NF}.lendis2| awk -v gt=$t -v rank=$count '{OFS="\t"}{ print gt".uniqmap",$1,$2,$3,rank}' >> ${lendisFile}
		count=$(($count+1))		 
		done
		echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_paired_lendis2 ${lendisFile} ${OUTDIR8}" >>${paraFile}
	done
	
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 4 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA.paired.lendis2WOroo
STEP=$((STEP+1))


declare -a GROUPGT=("pago3cdmut_ox" "pago3wtmut_ox" "pago3cdwt_ox" "pago3cdmut_unox" "pago3wtmut_unox" "pago3cdwt_unox" "ago3cdmut_ox" "ago3wtmut_ox" "ago3cdwt_ox" "ago3cdmut_unox" "ago3wtmut_unox" "ago3cdwt_unox" "aubcdmut_ox" "aubwtmut_ox" "aubcdwt_ox" "aubcdmut_unox" "aubwtmut_unox" "aubcdwt_unox" "ago3mut_cor1_ox" "ago3mut_cor1_unox" "ago3CD_cor2_ox" "ago3CD_cor2_unox")

declare -a pago3cdmut_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a pago3wtmut_ox=("Phil.SRA.nosAgo3WTrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a pago3cdwt_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.nosAgo3WTrescue.ox.ovary.inserts")

declare -a pago3cdmut_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a pago3wtmut_unox=("Phil.SRA.nosAgo3WTrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a pago3cdwt_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.nosAgo3WTrescue.unox.ovary.inserts")

declare -a ago3cdmut_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a ago3wtmut_ox=("Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a ago3cdwt_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts")

declare -a ago3cdmut_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a ago3wtmut_unox=("Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a ago3cdwt_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts")

declare -a aubcdmut_ox=("Phil.SRA.AubCDrescue.ox.ovary.inserts" "Phil.SRA.AubMutsWW.ox.ovary.inserts")
declare -a aubwtmut_ox=("Phil.SRA.AubWTrescue.ox.ovary.inserts" "Phil.SRA.AubMutsWW.ox.ovary.inserts")
declare -a aubcdwt_ox=("Phil.SRA.AubCDrescue.ox.ovary.inserts" "Phil.SRA.AubWTrescue.ox.ovary.inserts")

declare -a aubcdmut_unox=("Phil.SRA.AubCDrescue.unox.ovary.inserts" "Phil.SRA.AubMutsWW.unox.ovary.inserts")
declare -a aubwtmut_unox=("Phil.SRA.AubWTrescue.unox.ovary.inserts" "Phil.SRA.AubMutsWW.unox.ovary.inserts")
declare -a aubcdwt_unox=("Phil.SRA.AubCDrescue.unox.ovary.inserts" "Phil.SRA.AubWTrescue.unox.ovary.inserts")



declare -a ago3mut_cor1_ox=("Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.ago3MutsCJ.ox.ovary.inserts")
declare -a ago3mut_cor1_unox=("Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.ago3MutsCJ.unox.ovary.inserts")

declare -a ago3CD_cor2_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts")
declare -a ago3CD_cor2_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts")

echo -e "`date` "+$ISO_8601"\tDraw paired abundance,sense_fraction of transposon piRNAs" >> $LOG
OUTDIR9=${INDIR}/transposon_piRNA/paired_abundance_senseFraction
[ ! -d $OUTDIR9 ] && mkdir -p ${OUTDIR9}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.paired.abundance_senseFraction ] && \
paraFile=${OUTDIR9}/${RANDOM}.drawpairedabundance_senseFraction.para && \
for g in "${GROUPGT[@]}"
do
	SUBGROUP="$g[@]"
	for NF in "${NORMFACTORTYPE[@]}"
	do

		count=1
		transposonListFile=${OUTDIR9}/${g}.${NF}.transposon.list
		transposonListUniqFile=${OUTDIR9}/${g}.${NF}.uniqmap.transposon.list

		nnc=4 #the column in the stat file
		seqDep=2 #the column in the stat file
		for t in ${!SUBGROUP}
		do
			normFactor=`cat ${INDIR}/${t}/output/${t}_stats_table_reads|tail -1|cut -f${!NF}`
	
			tu=${t}.uniqmap
			normFactorUniq=`cat ${INDIR}/${t}/output/${tu}_stats_table_reads|tail -1|cut -f${!NF}`
	
			cat ${INDIR}/${t}/output/${t}.transposon.list| awk -v gt=$t -v rank=$count -v nf=$normFactor '{OFS="\t"}{if(rank==1) {tmp=match($0,/transposon/); if(tmp){print "gt",$0,"rank","nf"}else{print gt,$0,rank,nf}}else{tmp=match($0,/transposon/); if(!tmp){print gt,$0,rank,nf}}}' >> ${transposonListFile}
			cat ${INDIR}/${t}/output/${tu}.transposon.list| awk -v gt=$t -v rank=$count -v nf=$normFactorUniq '{OFS="\t"}{if(rank==1) { tmp=match($0,/transposon/); if(tmp) {print "gt",$0,"rank","nf"} else {print gt,$0,rank,nf} } else{tmp=match($0,/transposon/); if(!tmp){print gt,$0,rank,nf}}}' >> ${transposonListUniqFile}

			count=$(($count+1))
		done
		echo -e " ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_transposon_abundance_senseFraction_comparison ${transposonListFile} $NF ${OUTDIR9} " >>${paraFile}
		echo -e " ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_transposon_abundance_senseFraction_comparison ${transposonListUniqFile} $NF ${OUTDIR9} " >>${paraFile}
		
	done
	
done
[ $? == 0 ] && \
		ParaFly -c $paraFile -CPU 4 -failed_cmds $paraFile.failed_commands &&
	touch ${OUT}/.status.${STEP}.transposon_piRNA.paired.abundance_senseFraction
STEP=$((STEP+1))