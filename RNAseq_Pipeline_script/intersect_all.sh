#! /bin/bash
case ${3} in
## mouse specific intersectBed target files
"mouse")
	FOLDER=$DEG_PIPELINE_ADD/mm9/UCSC_annotation
;;
## fly specific intersectBed target files
"fly")
	FOLDER=$DEG_PIPELINE_ADD/dm3/UCSC_annotation
	
	PIRNA_CLUSTER=$FOLDER/Brennecke.pirnaCluster.bed
	REPEATMASKER_ALL=$FOLDER/repeat_mask.bed
	REPEATMASKER_DNA=$FOLDER/repeat_mask.bed.DNA
	REPEATMASKER_LINE=$FOLDER/repeat_mask.bed.LINE
	REPEATMASKER_LOW_COMPLEX=$FOLDER/repeat_mask.bed.Low_complexity
	REPEATMASKER_LTR=$FOLDER/repeat_mask.bed.LTR
	REPEATMASKER_OTHER=$FOLDER/repeat_mask.bed.Other
	REPEATMASKER_RC=$FOLDER/repeat_mask.bed.RC
	REPEATMASKER_RNA=$FOLDER/repeat_mask.bed.RNA
	REPEATMASKER_rRNA=$FOLDER/repeat_mask.bed.rRNA
	REPEATMASKER_SATELLITE=$FOLDER/repeat_mask.bed.Satellite
	REPEATMASKER_SIMPLE_REP=$FOLDER/repeat_mask.bed.Simple_repeat
	REPEATMASKER_UNKNOWN=$FOLDER/repeat_mask.bed.Unknown
	FLYBASE_3UTR=$FOLDER/UCSC.flyBase.3UTR.bed
	FLYBASE_5UTR=$FOLDER/UCSC.flyBase.5UTR.bed
	FLYBASE_EXON=$FOLDER/UCSC.flyBase.Exon.bed
	FLYBASE_GENE=$FOLDER/UCSC.flyBase.Genes.bed
	FLYBASE_INTRON=$FOLDER/UCSC.flyBase.Intron.bed

	declare -a TARGETS=("PIRNA_CLUSTER" \
     "REPEATMASKER_ALL" \
     "REPEATMASKER_DNA" \
     "REPEATMASKER_LINE" \
     "REPEATMASKER_LOW_COMPLEX" \
     "REPEATMASKER_LTR" \
     "REPEATMASKER_OTHER" \
     "REPEATMASKER_RC" \
     "REPEATMASKER_RNA" \
     "REPEATMASKER_rRNA" \
     "REPEATMASKER_SATELLITE" \
     "REPEATMASKER_SIMPLE_REP" \
     "REPEATMASKER_UNKNOWN" \
	"FLYBASE_3UTR" \
	"FLYBASE_5UTR" \
	"FLYBASE_EXON" \
    "FLYBASE_GENE" \
	"FLYBASE_INTRON" )
;;
*)
	echo "unknown orgnanism... currently only mouse/fly is supported"
	exit 2
;;
esac

UNIQUE_BED=$1
MULTIP_BED=$2
GENOME_FA=$4
# @ running pairtobed 
parafly_file="pairtobed".para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	echo "bedtools pairtobed -a ${UNIQUE_BED} -b ${!t} -bedpe -type ospan -f 0.75 > ${UNIQUE_BED%bam}${t}.bedpe" >> $parafly_file
	echo "bedtools pairtobed -a ${MULTIP_BED} -b ${!t} -bedpe -type ospan -f 0.75 > ${MULTIP_BED%bam}${t}.bedpe" >> $parafly_file
done
[ ! -f ${parafly_file}.completed ] && \
	ParaFly -c $parafly_file -CPU $CPU

# @ running bedpe S/AS counting
parafly_file="countPE".para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	echo "awk 'BEGIN{OFS=\"\\t\"} {if (ARGIND==1) {ct[\$4]=0; ct_pos[\$4]=0; ct_neg[\$4]=0;} if (ARGIND==2) {ct[\$14]++; if (\$9==\$16){ct_pos[\$14]++} else {ct_neg[\$14]++}}}END{for (i in ct){ print i,ct[i],ct_pos[i],ct_neg[i] }}' ${!t} ${UNIQUE_BED%bam}${t}.bedpe  > ${UNIQUE_BED%bam}${t}.bedpe.ct " >> $parafly_file && \
	echo "awk 'BEGIN{OFS=\"\\t\"} {if (ARGIND==1) {ct[\$4]=0; ct_pos[\$4]=0; ct_neg[\$4]=0;} if (ARGIND==2) {ct[\$14]++; if (\$9==\$16){ct_pos[\$14]++} else {ct_neg[\$14]++}}}END{for (i in ct){ print i,ct[i],ct_pos[i],ct_neg[i] }}' ${!t} ${MULTIP_BED%bam}${t}.bedpe  > ${MULTIP_BED%bam}${t}.bedpe.ct " >> $parafly_file
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

# @ running weblogo drawing
parafly_file="weblogo".para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	echo "bed22weblogo.sh ${UNIQUE_BED%bam}${t}.bedpe $GENOME_FA" >> $parafly_file && \
	echo "bed22weblogo.sh ${MULTIP_BED%bam}${t}.bedpe $GENOME_FA" >> $parafly_file
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi
