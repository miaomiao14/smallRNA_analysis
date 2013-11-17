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
OUT=${2}/materTables
LOG=${OUT}/log

declare -a NORMFACTORTYPE=("nnc" "seqDep")
nnc=4 #the column in the stat file
seqDep=2 #the column in the stat file

STEP=1

#get normalized total transposon mappers, transposon piRNAs, transposon sense piRNAs and transposon antisense piRNAs,sense fraction
echo -e "`date` "+$ISO_8601"\tDraw paired abundance,sense_fraction of transposon piRNAs" >> $LOG
OUTDIR=${OUT}/transposon_piRNAs_abundance_senseFraction
[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.abundance_senseFraction ] && \
paraFile=${OUTDIR}/${RANDOM}.abundance_senseFraction.para && \

for i in `ls ${INDIR}/*.inserts/output/*.transposon.list`
do 

	FILE=${i##*/}
	insertsname=`basename $FILE .transposon.list`
	for NF in "${NORMFACTORTYPE[@]}"
	do
	normFactor=`cat ${INDIR}/${t}/output/${t}_stats_table_reads|tail -1|cut -f${!NF}`
	colNum=`awk '{print NF}' ${i} | sort -nu | tail -n 1`
	done
#polarHistogram for sense fraction
${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/polarHistogram_sense_fraction.r plot_PolarHisto_senseFraction $i ${OUTDIR1}
done


for g in "${GROUPGT[@]}"
do
	SUBGROUP="$g[@]"
	for NF in "${NORMFACTORTYPE[@]}"
	do

		count=1
		transposonListFile=${OUTDIR}/${g}.${NF}.transposon.list
		rm ${OUTDIR}/${g}.${NF}.transposon.list
		transposonListUniqFile=${OUTDIR}/${g}.${NF}.uniqmap.transposon.list
		rm ${OUTDIR}/${g}.${NF}.uniqmap.transposon.list

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
		echo -e " ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_transposon_abundance_senseFraction_comparison ${transposonListFile} $NF ${OUTDIR} " >>${paraFile}
		echo -e " ${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_transposon_abundance_senseFraction_comparison ${transposonListUniqFile} $NF ${OUTDIR} " >>${paraFile}
	done
	
done
[ $? == 0 ] && \
		ParaFly -c $paraFile -CPU 4 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.transposon_piRNA.paired.abundance_senseFraction
STEP=$((STEP+1))

# transpose the table
ruby -lane 'BEGIN{$arr=[]}; $arr.concat([$F]); END{$arr.transpose.each{ |a| puts a.join ("\t") } }' -F"\t" $OUT > ${OUT}.t  && \
mv ${OUT}.t ${OUT} 