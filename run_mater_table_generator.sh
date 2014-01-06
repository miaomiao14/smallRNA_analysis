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
	
	for i in `ls ${INDIR}/*.inserts/output/*.transposon.list|grep -v uniqmap`
	do 
		
		FILE=${i##*/}
		insertsname=`basename $FILE .transposon.list` #genotype
		samplename=${insertsname#Phil.SRA.*}
		sample=${samplename/.ovary.inserts/}
		for NF in "${NORMFACTORTYPE[@]}"
		do
		normFactor=`cat ${INDIR}/${insertsname}/output/${insertsname}_stats_table_reads|tail -1|cut -f${!NF}`
		#colNum=`awk '{print NF}' ${i} | sort -nu | tail -n 1`
		awk -v nf=$normFactor -v gt=${sample} '{OFS="\t"}{print gt,$1,$2*1000000/nf,$3*1000000/nf,$4*1000000/nf,$5,$6*1000000/nf,$7*1000000/nf,$8*1000000/nf,$9,$10*1000000/nf,$11*1000000/nf,$12*1000000/nf,$13}' ${i} >${i%.transposon.list}.${NF}.normalized.transposon.list
		done
	done
	#ParaFly -c $paraFile -CPU 24 -failed_cmds $paraFile.failed_commands && \
	for NF in "${NORMFACTORTYPE[@]}"
	do
		for i in `ls ${INDIR}/*.inserts/output/*.${NF}.normalized.transposon.list`
		do
			COUNT=3
			for j in "${TRNLISTCOLNAMES[@]}"
			do
				cut -f1,2,${COUNT} ${i}|tail -n +2 >> ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.raw.txt
				COUNT=$(($COUNT+1))
			done
		done
	
	done
	for NF in "${NORMFACTORTYPE[@]}"
	do
		for j in "${TRNLISTCOLNAMES[@]}"
		do
			${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.raw.txt ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.txt
		done
	done
	
	for NF in "${NORMFACTORTYPE[@]}"
	do
		for j in "${TRNLISTCOLNAMES[@]}"
		do
			${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source grouping ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.txt ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.withgroup.txt && rm ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.raw.txt
		done
	done
fi



[ $? == 0 ] && \
	#ParaFly -c $paraFile -CPU 24 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.transposon_piRNA.paired.abundance_senseFraction
STEP=$((STEP+1))
echo -e "`date` "+$ISO_8601"\tcreate master table for all genotypes abundance,sense_fraction of transposon piRNAs done" >> $LOG
# transpose the table
#ruby -lane 'BEGIN{$arr=[]}; $arr.concat([$F]); END{$arr.transpose.each{ |a| puts a.join ("\t") } }' -F"\t" $OUT > ${OUT}.t  && \
#mv ${OUT}.t ${OUT} 


#get normalized(to siRNAs map to structure loci and cisNats) total transposon mappers, transposon piRNAs, transposon sense piRNAs and transposon antisense piRNAs,sense fraction
GT=("w1.ox" "ago3Muts.ox" "ago3Mutsrep2.ox" "aubvasAgo3WTrescue.ox" "aubvasAgo3WTrescuerep2.ox" "aubvasAgo3CDrescue.ox" "aubvasAgo3CDrescuerep2.ox" "AubMutsWW.ox" "aubMutsrep2.ox" "AubWTrescue.ox" "AubWTrescuerep2.ox" "AubCDrescue.ox" "AubCDrescuerep2.ox" "AubCDinHets.ox" "ago3MutsAubMuts.ox" "QinHets.ox" "QinMuts.ox" "QinAgo3Hets.ox" "QinAgo3Muts.ox")
OUTDIR=${OUT}/nf_siRNA_transposon_piRNAs_abundance_senseFraction
NFFILE=/home/wangw1/pipeline_dm3/common/nf_siRNA.txt
NF="siRNA"
[ ! -d $OUTDIR ] && mkdir -p ${OUTDIR}
if [ ! -f ${OUT}/.status.${STEP}.nf_siRNA_transposon_piRNA.abundance_senseFraction ]
then
	for i in ${GT[@]}
	do
		file=${INDIR}/Phil.SRA.${i}.ovary.inserts/output/Phil.SRA.${i}.ovary.inserts.transposon.list
		FILE=${file##*/}
		insertsname=`basename $FILE .transposon.list` #genotype
		
		sample=$i
		normFactor=`grep $i $NFFILE|cut -f2`
		awk -v nf=$normFactor -v gt=${sample} '{OFS="\t"}{print gt,$1,$2/nf,$3/nf,$4/nf,$5,$6/nf,$7/nf,$8/nf,$9,$10/nf,$11/nf,$12/nf,$13}' ${file} >${file%.transposon.list}.${NF}.normalized.transposon.list
		
	done
	
	for i in ${GT[@]}
	do
		file=${INDIR}/Phil.SRA.${i}.ovary.inserts/output/Phil.SRA.${i}.ovary.inserts.transposon.list
		COUNT=3
		for j in "${TRNLISTCOLNAMES[@]}"
		do
			cut -f1,2,${COUNT} ${file%.transposon.list}.${NF}.normalized.transposon.list|tail -n +2 >> ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.raw.txt
			COUNT=$(($COUNT+1))
		done
	done
	for j in "${TRNLISTCOLNAMES[@]}"
	do
		${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source grouping ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.txt ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.withgroup.txt && rm ${OUTDIR}/${j}.${NF}.normalized.SRA.TRN.mastertable.raw.txt
	done
	
fi
[ $? == 0 ] && \
	#ParaFly -c $paraFile -CPU 24 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.nf_siRNA_transposon_piRNA.paired.abundance_senseFraction
STEP=$((STEP+1))
echo -e "`date` "+$ISO_8601"\tcreate master table for all ox genotypes abundance,sense_fraction of transposon piRNAs normalize to siRNAs done" >> $LOG


#master table for DEG htseq-count(including multiple mappers)

#export RUN_PIPELINE_ADD=/home/wangw1/git/PE_RNA_Pipeline
#Genome="dm3"
#declare -a TARGETS=("FLY_PIRNA_CLUSTER" \
#	"FLY_TRANSPOSON_ALL" \
#	"FLY_GENE_TRANSPOSON_ALL" \
#	"FLY_TRANSPOSON_OUTCLUSTER" \
#	"FLY_TRANSPOSON_INCLUSTER" \ 
#	"FLY_TRANSPOSON_GROUP1" \
#	"FLY_TRANSPOSON_GROUP2" \
#	"FLY_TRANSPOSON_GROUP3" \
#	"FLY_TRANSPOSON_GROUP0" \
#     "FLY_flyBase_GENE" \
#     "FLY_flyBase_EXON" \
#     "FLY_flyBase_INTRON" \
#     "FLY_flyBase_5UTR" \
#     "FLY_flyBase_3UTR" )
#    
#    
#declare -a GROUPGT=("GROUPGT_R1")
#
#declare -a GROUPGT_R1=("ago3cdmut_r1" "ago3wtmut_r1" "ago3cdwt_r1" "aubcdmut_r1" "aubwtmut_r1" "aubcdwt_r1")
#declare -a ago3cdmut_r1=("aubvasAgo3CDrescue" "ago3MutsWW")
#declare -a ago3wtmut_r1=("aubvasAgo3WTrescue" "ago3MutsWW")
#declare -a ago3cdwt_r1=("aubvasAgo3CDrescue" "aubvasAgo3WTrescue")
#declare -a aubcdmut_r1=("AubCDrescue" "aubMutsWW")
#declare -a aubwtmut_r1=("AubWTrescue" "aubMutsWW")
#declare -a aubcdwt_r1=("AubCDrescue" "AubWTrescue")
#

