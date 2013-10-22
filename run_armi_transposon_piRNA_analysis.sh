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
BINSIZE=$2 #BINSIZE for cicos plot
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log
#
declare -a GROUPGT=("wtk729a_ox" "wthets_ox" "wtmut_ox" "k729amut_ox" "k729ahets_ox" "k729aPNK_unox" "acthets_ox" "muthets_ox" "copy1vscopy2_ox")
declare -a wtk729a_ox=("Phil.SRA.nosactArmiWTrescue72D2.ox.ovary.inserts" "Phil.SRA.nosactArmiK729Arescue72D2.ox.ovary.inserts")
declare -a wthets_ox=("Phil.SRA.nosactArmiWTrescue72D2.ox.ovary.inserts" "Phil.SRA.nosactinarmiD2Hets.ox.ovary.inserts")
declare -a wtmut_ox=("Phil.SRA.nosactArmiWTrescue72D2.ox.ovary.inserts" "Phil.SRA.ArmiWTin72D2.ox.ovary.inserts")
declare -a k729amut_ox=("Phil.SRA.nosactArmiK729Arescue72D2.ox.ovary.inserts" "Phil.SRA.ArmiWTin72D2.ox.ovary.inserts")
declare -a k729ahets_ox=("Phil.SRA.nosactArmiK729Arescue72D2.ox.ovary.inserts" "Phil.SRA.nosactinarmiD2Hets.ox.ovary.inserts")
declare -a k729aPNK_unox=("Phil.SRA.nosactArmiK729Arescue72D2.unox.ovary.inserts" "Phil.SRA.nosactArmiK729Arescue72D2_PNK.unox.ovary.inserts")
declare -a acthets_ox=("Phil.SRA.armi72_1Hets.ox.ovary.inserts" "Phil.SRA.nosactinarmiD2Hets.ox.ovary.inserts")
declare -a muthets_ox=("Phil.SRA.ArmiWTin72D2.ox.ovary.inserts" "Phil.SRA.nosactinarmiD2Hets.ox.ovary.inserts")
declare -a copy1vscopy2_ox=("Phil.SRA.nosactinarmiD2Hets.ox.ovary.inserts" "Phil.SRA.UASPArmiWTinHets.ox.ovary.inserts")


#declare -a GROUPGT=("ago3hetmut_piwiip_unox")

#declare -a ago3hetmut_piwiip_unox=("Phil.SRA.PiwiIPago3hets.unox.ovary.inserts" "Phil.SRA.PiwiIPago3muts.unox.ovary.inserts")

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
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.transposon_piRNA.lendis2
STEP=$((STEP+1))



#declare -a gt_g1_ox=("aubvasAgo3CDrescue.ox" "ago3MutsWW.ox" "aubvasAgo3WTrescue.ox")
#declare -a gt_g1_unox=("aubvasAgo3CDrescue.unox" "ago3MutsWW.unox" "aubvasAgo3WTrescue.unox")
#declare -a gt_g2_ox=("AubCDrescue.ox" "aubMutsWW.ox" "AubWTrescue.ox")
#declare -a gt_g2_unox=("AubCDrescue.unox" "aubMutsWW.unox" "AubWTrescue.unox")
#
#declare -a gt_cor1_ox=("ago3MutsWW.ox" "ago3MutsCJ.ox")
#declare -a gt_cor1_unox=("ago3MutsWW.unox" "ago3MutsCJ.unox")
#
#declare -a gt_cor2_ox=("nosAgo3CDrescue.ox" "aubvasAgo3CDrescue.ox")
#declare -a gt_cor2_unox=("nosAgo3CDrescue.unox" "aubvasAgo3CDrescue.unox")

#declare -a GROUPGT=("armi_g1_ox" "armi_g1_unox" \
#"armi_g2_ox" "armi_g2_unox" \
#"armi_g3_ox" "armi_g3_unox" \
#"armi_wtk729a_ox" "armi_wtk729a_unox" \
#"armi_wtrescue_cor1" "armi_k729arescue_cor1" \
#"armi_armi72_1Hets_cor1" "armi_wtinhets_cor1" \
#)
#
#declare -a armi_g1_ox=("Phil.SRA.nosGal4UASPArmiWTrescue.ox.ovary.inserts" "Phil.SRA.UASPArmiWTinHets.ox.ovary.inserts")
#declare -a armi_g1_unox=("Phil.SRA.nosGal4UASPArmiWTrescue.unox.ovary.inserts" "Phil.SRA.UASPArmiWTinHets.unox.ovary.inserts")
#
#declare -a armi_g2_ox=("Phil.SRA.nosGal4UASPArmiK729Arescue.ox.ovary.inserts" "Phil.SRA.UASPArmiWTinHets.ox.ovary.inserts")
#declare -a armi_g2_unox=("Phil.SRA.nosGal4UASPArmiK729Arescue.unox.ovary.inserts" "Phil.SRA.UASPArmiWTinHets.unox.ovary.inserts")
#
#declare -a armi_g3_ox=("Phil.SRA.armi72_1Hets.ox.ovary.inserts" "Phil.SRA.UASPArmiWTinHets.ox.ovary.inserts")
#declare -a armi_g3_unox=("Phil.SRA.armi72_1Hets.unox.ovary.inserts" "Phil.SRA.UASPArmiWTinHets.unox.ovary.inserts")
#
#declare -a armi_wtk729a_ox=("Phil.SRA.nosGal4UASPArmiWTrescue.ox.ovary.inserts" "Phil.SRA.nosGal4UASPArmiK729Arescue.ox.ovary.inserts")
#declare -a armi_wtk729a_unox=("Phil.SRA.nosGal4UASPArmiWTrescue.unox.ovary.inserts" "Phil.SRA.nosGal4UASPArmiK729Arescue.unox.ovary.inserts")
#
#declare -a armi_wtrescue_cor1=("Phil.SRA.nosGal4UASPArmiWTrescue.ox.ovary.inserts" "Phil.SRA.nosGal4UASPArmiWTrescue.unox.ovary.inserts")
#declare -a armi_k729arescue_cor1=("Phil.SRA.nosGal4UASPArmiK729Arescue.ox.ovary.inserts" "Phil.SRA.nosGal4UASPArmiK729Arescue.unox.ovary.inserts")
#declare -a armi_armi72_1Hets_cor1=("Phil.SRA.armi72_1Hets.ox.ovary.inserts" "Phil.SRA.armi72_1Hets.unox.ovary.inserts")
#declare -a armi_wtinhets_cor1=("Phil.SRA.UASPArmiWTinHets.ox.ovary.inserts" "Phil.SRA.UASPArmiWTinHets.unox.ovary.inserts")



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
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands && \
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
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands && \
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






INDIR6=${INDIR}/circos
OUTDIR6=/home/wangw1/src/circos-0.56/fly
OUTDIR6_2=${INDIR}/transposon_piRNA/circos_conf
#/home/wangw1/src/circos-0.56/bin/circos -conf /home/wangw1/src/circos-0.56/fly/etc/

declare -a NORMFACTORTYPE=("nnc" "seqDep")
declare -a METRICTYPE=("max" "mean")
[ ! -d $OUTDIR6/armicircosPlot ] && mkdir -p ${OUTDIR6}/armicircosPlot
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
#			echo -ne " mv ${OUTDIR6}/circos.svg ${OUTDIR6}/armicircosPlot/${g}.${NF}.${MT}.${BINSIZE}.svg && " >>${paraFile}
#			echo -e " mv ${OUTDIR6}/circos.png ${OUTDIR6}/armicircosPlot/${g}.${NF}.${MT}.${BINSIZE}.png " >>${paraFile}
#			
#			echo -ne " ${PIPELINE_DIRECTORY}/circos_conf_gnr.sh ${TARGETSUNIQ[@]} ${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			#echo -ne " mv ${OUTDIR6}/etc/file.conf ${OUTDIR6}/etc/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			echo -ne " cd ${OUTDIR6} && " >>${paraFile}
#			echo -ne " /home/wangw1/src/circos-0.56/bin/circos -conf ./etc/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
#			echo -ne " mv ${OUTDIR6}/circos.svg ${OUTDIR6}/armicircosPlot/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.svg && " >>${paraFile}
#			echo -e " mv ${OUTDIR6}/circos.png ${OUTDIR6}/armicircosPlot/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.png " >>${paraFile}
			
			${PIPELINE_DIRECTORY}/circos_conf_gnr.sh ${TARGETS[@]} ${g}.${NF}.${MT}.${BINSIZE}.conf && \

			cd ${OUTDIR6} && \
			/home/wangw1/src/circos-0.56/bin/circos -conf ./etc/${g}.${NF}.${MT}.${BINSIZE}.conf && \
			mv ${OUTDIR6}/circos.svg ${OUTDIR6}/armicircosPlot/${g}.${NF}.${MT}.${BINSIZE}.svg && \
			mv ${OUTDIR6}/circos.png ${OUTDIR6}/armicircosPlot/${g}.${NF}.${MT}.${BINSIZE}.png
			
			${PIPELINE_DIRECTORY}/circos_conf_gnr.sh ${TARGETSUNIQ[@]} ${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && \
			#echo -ne " mv ${OUTDIR6}/etc/file.conf ${OUTDIR6}/etc/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && " >>${paraFile}
			cd ${OUTDIR6} && \
			/home/wangw1/src/circos-0.56/bin/circos -conf ./etc/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.conf && \
			mv ${OUTDIR6}/circos.svg ${OUTDIR6}/armicircosPlot/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.svg && \
			mv ${OUTDIR6}/circos.png ${OUTDIR6}/armicircosPlot/${g}.uniqmap.${NF}.${MT}.${BINSIZE}.png
			
			
		done
	done
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.transposon_piRNA.circosplot
STEP=$((STEP+1))



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
		rm ${OUTDIR9}/${g}.${NF}.transposon.list
		transposonListUniqFile=${OUTDIR9}/${g}.${NF}.uniqmap.transposon.list
		rm ${OUTDIR9}/${g}.${NF}.uniqmap.transposon.list

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
		ParaFly -c $paraFile -CPU 4 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.transposon_piRNA.paired.abundance_senseFraction
STEP=$((STEP+1))



#touch ${OUT}/.status.${STEP}.transposon_piRNA.pairedZscore
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

#transposon bucket
echo -e "`date` "+$ISO_8601"\trun bucket of transposon piRNAs" >> $LOG
OUTDIR11=${INDIR}/transposon_piRNA/transposon_bucket
[ ! -d $OUTDIR11 ] && mkdir -p ${OUTDIR11}
[ ! -f ${OUT}/.status.${STEP}.transposon_bucket ] && \
for g in "${GROUPGT[@]}"
do
	#SUBGROUP="$g[@]"
		eval "SUBGROUP=(\"\${${g}[@]}\")"  #array in bash can not be assigned directly

		[ ! -d ${OUTDIR11}/${g} ] && mkdir ${OUTDIR11}/${g}
		outputdir=${OUTDIR11}/${g}
		inputdir=$outputdir
		inputfilename1=${SUBGROUP[0]}
		inputfilename2=${SUBGROUP[1]}
		ln -s ${INDIR}/${inputfilename1}/${inputfilename1}.xkxh.transposon.mapper2.gz ${outputdir}
		ln -s ${INDIR}/${inputfilename2}/${inputfilename2}.xkxh.transposon.mapper2.gz ${outputdir}
		samplename1b=${inputfilename1#Phil.SRA.*}
		samplename2b=${inputfilename2#Phil.SRA.*}
		samplename1=${samplename1b%*.ox.ovary.inserts}
		samplename2=${samplename2b%*.ox.ovary.inserts}
		seqdepth1=`cat ${INDIR}/${inputfilename1}/output/${inputfilename1}_stats_table_reads|tail -1|awk '{print $4/1000000}'`
		seqdepth2=`cat ${INDIR}/${inputfilename2}/output/${inputfilename2}_stats_table_reads|tail -1|awk '{print $4/1000000}'`
		email="weiwanghhq@gmail.com"
		
${PIPELINE_DIRECTORY}/bucket_new_gz_batch.pl ${inputfilename1}.xkxh.transposon.mapper2.gz ${inputfilename2}.xkxh.transposon.mapper2.gz ${inputdir} ${outputdir} ${samplename1} ${samplename2} ${seqdepth1} ${seqdepth2} ${email} >$LOG
	
	
done
touch ${OUT}/.status.${STEP}.transposon_bucket
STEP=$((STEP+1))

#cluster bucket
echo -e "`date` "+$ISO_8601"\trun cluster bucket of transposon piRNAs" >> $LOG
OUTDIR12=${INDIR}/transposon_piRNA/cluster_bucket
[ ! -d $OUTDIR12 ] && mkdir -p ${OUTDIR12}
# USAGE: make_sge_cluster_bucket.sh <n> <input xkxh.norm.bed file1 with full path>  <stats table1> <input xkxh.norm.bed file2 with full path> <stats table2> <n> <option [empty: default Brennecke 142 ovary clusters;  custom defined cluster file]>
[ ! -f ${OUT}/.status.${STEP}.cluster_bucket ] && \
for g in "${GROUPGT[@]}"
do
	eval "SUBGROUP=(\"\${${g}[@]}\")"
	[ ! -d ${OUTDIR12}/${g} ] && mkdir ${OUTDIR12}/${g}
	outputdir=${OUTDIR12}/${g}
	inputdir=$outputdir
	inputfilename1=${SUBGROUP[0]}
	inputfilename2=${SUBGROUP[1]}
	#ln -s ${INDIR}/${inputfilename1}/${inputfilename1}.xkxh.norm.bed.gz ${outputdir}
	#ln -s ${INDIR}/${inputfilename2}/${inputfilename2}.xkxh.norm.bed.gz ${outputdir}
	stattable1=${INDIR}/${inputfilename1}/output/${inputfilename1}_stats_table_reads
	stattable2=${INDIR}/${inputfilename2}/output/${inputfilename2}_stats_table_reads
	${PIPELINE_DIRECTORY}/submit_cluster_bucket_ww.sh 2 ${INDIR}/${inputfilename1}/${inputfilename1}.xkxh.norm.bed.gz ${stattable1} ${INDIR}/${inputfilename2}/${inputfilename2}.xkxh.norm.bed.gz ${stattable2}
done
touch ${OUT}/.status.${STEP}.cluster_bucket


echo -e "`date` "+$ISO_8601"\tDraw phasing analysis..." >> $LOG
OUTDIR13=${INDIR}/transposon_piRNA/phasing
[ ! -d $OUTDIR13 ] && mkdir -p ${OUTDIR13}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.phasing ] && \
paraFile=${OUTDIR13}/${RANDOM}.piRNAphasing.para

#for i in `ls ${INDIR}/*.inserts`
#do
	#	inputfile=
#${PIPELINE_DIRECTORY}/run_distance_analysis.sh
