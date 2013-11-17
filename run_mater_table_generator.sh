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
OUT=${2}/masterTables
LOG=${OUT}/log

declare -a NORMFACTORTYPE=("nnc")
nnc=4 #the column in the stat file
seqDep=2 #the column in the stat file

STEP=1

#get normalized total transposon mappers, transposon piRNAs, transposon sense piRNAs and transposon antisense piRNAs,sense fraction

TRNLISTCOLNAMES=("total" "total_sense" "total_antisense" "total_sense_fraction" "siRNA" "siRNA_sense" "siRNA_antisense" "siRNA_sense_fraction" "piRNA" "piRNA_sense" "piRNA_antisense" "piRNA_sense_fraction")
echo -e "`date` "+$ISO_8601"\tcreate master table for all genotypes abundance,sense_fraction of transposon piRNAs" >> $LOG
OUTDIR=${OUT}/transposon_piRNAs_abundance_senseFraction
[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}
if [ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.abundance_senseFraction ]
then
	paraFile=${OUTDIR}/${RANDOM}.abundance_senseFraction.para && \
	
	for i in `ls ${INDIR}/*.inserts/output/*.transposon.list `
	do 
		
		FILE=${i##*/}
		insertsname=`basename $FILE .transposon.list` #genotype
		for NF in "${NORMFACTORTYPE[@]}"
		do
		normFactor=`cat ${INDIR}/${insertsname}/output/${insertsname}_stats_table_reads|tail -1|cut -f${!NF}`
		#colNum=`awk '{print NF}' ${i} | sort -nu | tail -n 1`
		awk -v nf=1000000/$normFactor -v gt=${insertsname} '{OFS="\t"}{print gt,$1,$2*nf,$3*nf,$4*nf,$5,$6*nf,$7*nf,$8*nf,$9,$10*nf,$11*nf,$12*nf,$13}' ${i} >${i%.transposon.list}.${NF}.normalized.transposon.list
		done
	done
	ParaFly -c $paraFile -CPU 24 -failed_cmds $paraFile.failed_commands && \
	for NF in "${NORMFACTORTYPE[@]}"
	do
		for i in `ls ${INDIR}/*.inserts/output/*.${NF}.normalized.transposon.list`
		do
			COUNT=3
			for j in "${TRNLISTCOLNAMES[@]}"
			do
				cut -f1,2,${COUNT} ${i}|tail -n +2 >> ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.txt
				COUNT=$(($COUNT+1))
			done
		done
	
	done
fi
#${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/polarHistogram_sense_fraction.r plot_PolarHisto_senseFraction $i ${OUTDIR1}

[ $? == 0 ] && \
	#ParaFly -c $paraFile -CPU 24 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.transposon_piRNA.paired.abundance_senseFraction
STEP=$((STEP+1))
echo -e "`date` "+$ISO_8601"\tcreate master table for all genotypes abundance,sense_fraction of transposon piRNAs done" >> $LOG
# transpose the table
#ruby -lane 'BEGIN{$arr=[]}; $arr.concat([$F]); END{$arr.transpose.each{ |a| puts a.join ("\t") } }' -F"\t" $OUT > ${OUT}.t  && \
#mv ${OUT}.t ${OUT} 