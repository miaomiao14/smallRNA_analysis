#! /bin/bash -x

#RUN DIRECTORY
export RUN_PIPELINE_ADD=/home/wangw1/git/PE_RNA_Pipeline

case ${3} in
## mouse specific intersectBed target files
"mouse")
	FOLDER=/home/hanb/nearline/small_RNA_Pipeline/common_files/mm9/UCSC_BEDS
	# piRNA clusters annoated in Li, et al, Mol Cell, 2013
	MOUSE_PIRNACLUSTER=$FOLDER/piRNA.cluster.bed
	MOUSE_PREPACHYTENE_PIRNA_CLUSTER=$FOLDER/piRNA.cluster.prepachytene.bed
	MOUSE_HYBRID_PIRNA_CLUSTER=$FOLDER/piRNA.cluster.hybrid.bed
	MOUSE_PACHYTENE_PIRNA_CLUSTER=$FOLDER/piRNA.cluster.pachytene.bed
	# NM transcripts annotated in Zamore lab, 2013
	MOUSE_ZAMORE_NM=$FOLDER/x100910.rnaseq.transcripts.NM.all.final.bed.c
	# NR transcripts annotated in Zamore lab, 2013
	MOUSE_ZAMORE_NR=$FOLDER/x100910.rnaseq.transcripts.NR.all.final.bed.c
	# transcriptome combined from piRNA cluster, NM and NR
	MOUSE_ZAMORE_TRANSCRIPTOME=$FOLDER/zamore.transcriptome.bed
	# take the col1-6 of $MOUSE_ZAMORE_TRANSCRIPTOME (including intron)
	MOUSE_ZAMORE_ANNOTATION=$FOLDER/zamore.annotation.bed 
	
	MOUSE_refSeq_GENE=$FOLDER/UCSC.refSeq.Genes.bed
	MOUSE_refSeq_EXON=$FOLDER/UCSC.refSeq.Exon.bed
	MOUSE_refSeq_INTRON=$FOLDER/UCSC.refSeq.Intron.bed
	MOUSE_refSeq_5UTR=$FOLDER/UCSC.refSeq.5UTR.bed
	MOUSE_refSeq_3UTR=$FOLDER/UCSC.refSeq.3UTR.bed
	
	MOUSE_REPEATMASKER_DNA=$FOLDER/repeat_mask.bed.DNA
	MOUSE_REPEATMASKER_SINE=$FOLDER/repeat_mask.bed.SINE
	MOUSE_REPEATMASKER_LINE=$FOLDER/repeat_mask.bed.LINE
	MOUSE_REPEATMASKER_LTR=$FOLDER/repeat_mask.bed.LTR
	MOUSE_REPEATMASKER_SATELLITE=$FOLDER/repeat_mask.bed.Satellite
	MOUSE_REPEATMASKER_SIMPLEREPEATS=$FOLDER/repeat_mask.bed.Simple_repeat
	MOUSE_REPEATMASKER_tRNA=$FOLDER/repeat_mask.bed.tRNA
	
	declare -a TARGETS=( \
	"MOUSE_PIRNACLUSTER" \
    "MOUSE_PREPACHYTENE_PIRNA_CLUSTER" \
    "MOUSE_HYBRID_PIRNA_CLUSTER" \
    "MOUSE_PACHYTENE_PIRNA_CLUSTER" \
    "MOUSE_ZAMORE_NM" \
    "MOUSE_ZAMORE_NR" \
    "MOUSE_ZAMORE_TRANSCRIPTOME" \
    "MOUSE_refSeq_GENE" \
    "MOUSE_refSeq_EXON" \
    "MOUSE_refSeq_INTRON" \
    "MOUSE_refSeq_5UTR" \
    "MOUSE_refSeq_3UTR" \
    "MOUSE_REPEATMASKER_DNA" \
    "MOUSE_REPEATMASKER_SINE" \
    "MOUSE_REPEATMASKER_LINE" \
    "MOUSE_REPEATMASKER_LTR" \
    "MOUSE_REPEATMASKER_SATELLITE" \
    "MOUSE_REPEATMASKER_SIMPLEREPEATS" \
    "MOUSE_REPEATMASKER_tRNA" )
;;
## fly specific intersectBed target files
"fly")
	FOLDER=/home/hanb/nearline/small_RNA_Pipeline/common_files/dm3/UCSC_BEDS
	
	FLY_PIRNA_CLUSTER=$FOLDER/Brennecke.pirnaCluster.bed 
	FLY_PIRNA_CLUSTER_42AB=$FOLDER/Brennecke.pirnaCluster.42AB.bed
	FLY_PIRNA_CLUSTER_FLAM=$FOLDER/Brennecke.pirnaCluster.flam.bed
	FLY_flyBase_GENE=$FOLDER/UCSC.flyBase.Genes.bed
	FLY_flyBase_EXON=$FOLDER/UCSC.flyBase.Exon.bed
	FLY_flyBase_INTRON=$FOLDER/UCSC.flyBase.Intron.bed
	FLY_flyBase_INTRON_xRM=$FOLDER/UCSC.flyBase.Intron.xRM.bed
	FLY_flyBase_5UTR=$FOLDER/UCSC.flyBase.5UTR.bed
	FLY_flyBase_3UTR=$FOLDER/UCSC.flyBase.3UTR.bed
	FLY_REPEATMASKER=$FOLDER/repeat_mask.bed
	FLY_REPEATMASKER_IN_CLUSTER=$FOLDER/repeat_mask.inCluster.bed
	FLY_REPEATMASKER_OUT_CLUSTER=$FOLDER/repeat_mask.outCluster.bed
	FLY_REPEATMASKER_DNA=$FOLDER/repeat_mask.bed.DNA
	FLY_REPEATMASKER_LINE=$FOLDER/repeat_mask.bed.LINE
	FLY_REPEATMASKER_LTR=$FOLDER/repeat_mask.bed.LTR
	FLY_REPEATMASKER_RNA=$FOLDER/repeat_mask.bed.RNA
	FLY_REPEATMASKER_rRNA=$FOLDER/repeat_mask.bed.rRNA
	FLY_REPEATMASKER_Satellite=$FOLDER/repeat_mask.bed.Satellite
	FLY_REPEATMASKER_Simple_repeat=$FOLDER/repeat_mask.bed.Simple_repeat
	FLY_REPEATMASKER_Unknown=$FOLDER/repeat_mask.bed.Unknown
	
	FLY_TRANSPOSON_ALL=$FOLDER/transposon.bed2
	FLY_TRANSPOSON_GROUP1=$FOLDER/Zamore.group1.bed
	FLY_TRANSPOSON_GROUP2=$FOLDER/Zamore.group2.bed
	FLY_TRANSPOSON_GROUP3=$FOLDER/Zamore.group3.bed
	FLY_TRANSPOSON_GROUP0=$FOLDER/Zamore.group0.bed
	
	declare -a TARGETS=("FLY_PIRNA_CLUSTER" \
	"FLY_PIRNA_CLUSTER_42AB" \
	"FLY_PIRNA_CLUSTER_FLAM" \
	"FLY_TRANSPOSON_ALL" \
	"FLY_TRANSPOSON_GROUP1" \
	"FLY_TRANSPOSON_GROUP2" \
	"FLY_TRANSPOSON_GROUP3" \
	"FLY_TRANSPOSON_GROUP0" \
    "FLY_REPEATMASKER" \
    "FLY_REPEATMASKER_IN_CLUSTER" \
    "FLY_REPEATMASKER_OUT_CLUSTER" \
     "FLY_REPEATMASKER_DNA" \
     "FLY_REPEATMASKER_LINE" \
    "FLY_REPEATMASKER_LTR" \
     "FLY_flyBase_GENE" \
     "FLY_flyBase_EXON" \
     "FLY_flyBase_INTRON" \
    "FLY_flyBase_INTRON_xRM" \
     "FLY_flyBase_5UTR" \
     "FLY_flyBase_3UTR" \
    )
    
    FOLDER=/home/wangw1/pipeline_dm3/common
    FLY_PIRNA_CLUSTER=$FOLDER/Brennecke.pirnaCluster.bed6 
	FLY_flyBase_GENE=$FOLDER/UCSC.flyBase.Genes.bed6
	FLY_flyBase_EXON=$FOLDER/UCSC.flyBase.Exon.bed6
	FLY_flyBase_INTRON=$FOLDER/UCSC.flyBase.Intron.bed6
	FLY_flyBase_5UTR=$FOLDER/UCSC.flyBase.5UTR.bed6
	FLY_flyBase_3UTR=$FOLDER/UCSC.flyBase.3UTR.bed6
	FLY_REPEATMASKER_DNA=$FOLDER/repeat_mask.bed.DNA
	FLY_REPEATMASKER_LINE=$FOLDER/repeat_mask.bed.LINE
	FLY_REPEATMASKER_LTR=$FOLDER/repeat_mask.bed.LTR
	FLY_REPEATMASKER_RNA=$FOLDER/repeat_mask.bed.RNA
	FLY_REPEATMASKER_rRNA=$FOLDER/repeat_mask.bed.rRNA
	FLY_REPEATMASKER_Satellite=$FOLDER/repeat_mask.bed.Satellite
	FLY_REPEATMASKER_Simple_repeat=$FOLDER/repeat_mask.bed.Simple_repeat
	FLY_REPEATMASKER_Unknown=$FOLDER/repeat_mask.bed.Unknown
	
	FLY_TRANSPOSON_ALL=$FOLDER/transposon.group.bed6
	FLY_TRANSPOSON_GROUP1=$FOLDER/Zamore.group1.bed6
	FLY_TRANSPOSON_GROUP2=$FOLDER/Zamore.group2.bed6
	FLY_TRANSPOSON_GROUP3=$FOLDER/Zamore.group3.bed6
	FLY_TRANSPOSON_GROUP0=$FOLDER/Zamore.group0.bed6
	
	declare -a TARGETSBED6=("FLY_PIRNA_CLUSTER" \
	"FLY_TRANSPOSON_ALL" \
	"FLY_TRANSPOSON_GROUP1" \
	"FLY_TRANSPOSON_GROUP2" \
	"FLY_TRANSPOSON_GROUP3" \
	"FLY_TRANSPOSON_GROUP0" \
     "FLY_flyBase_GENE" \
     "FLY_flyBase_EXON" \
     "FLY_flyBase_INTRON" \
     "FLY_flyBase_5UTR" \
     "FLY_flyBase_3UTR" \
     "FLY_REPEATMASKER_DNA" \
     "FLY_REPEATMASKER_LINE" \
     "FLY_REPEATMASKER_LTR" \
     "FLY_REPEATMASKER_Satellite" \
     "FLY_REPEATMASKER_Simple_repeat" \
     "FLY_REPEATMASKER_RNA" \
     "FLY_REPEATMASKER_rRNA" \
     "FLY_REPEATMASKER_Unknown" )
    
;;
*)
	echo "unknown orgnanism... currently only mouse/fly is supported"
	exit 2
;;
esac

UNIQUE_BED=$1
MULTIP_BED=$2
UNIQUE_BED_NAME=${1##*/}
MULTIP_BED_NAME=${2##*/}
GENOME_FA=$4


PREFIX=${UNIQUE_BED_NAME%.uniq*}

INTERSECTOUTDIR=$6
LOGODIR=$7
CPU=8
#${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bam}
#${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bam}

#collapse reads to species with read counts,make sure no soft clip of the 5'end
sort -k1,1 -k2,2n -k3,3n -k5,5n ${UNIQUE_BED} | bedtools groupby -i - -g 1,2,3,4,5,6,8,9,10 -c 7,8 -o collapse,sum |awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$10,$11,$8,$9}' >${UNIQUE_BED}.collapse
sort -k1,1 -k2,2n -k3,3n -k5,5n ${MULTIP_BED} | bedtools groupby -i - -g 1,2,3,4,5,6,8,9,10 -c 7,8 -o collapse,sum |awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$10,$11,$8,$9}' >${MULTIP_BED}.collapse

# @ running pairtobed 
parafly_file=${INTERSECTOUTDIR}/pairtobed.para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	echo "bedtools pairtobed -a ${UNIQUE_BED}.collapse -b ${!t} -bedpe -type ospan -f 0.75 > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bam}${t}.bedpe" >> $parafly_file
	echo "bedtools pairtobed -a ${MULTIP_BED}.collapse -b ${!t} -bedpe -type ospan -f 0.75 > ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bam}${t}.bedpe" >> $parafly_file
	#echo "bedtools pairtobed -a ${UNIQUE_BED} -b ${!t} -bedpe -type ospan -f 0.75 -s > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bam}${t}.bedpe" >> $parafly_file
	#echo "bedtools pairtobed -a ${MULTIP_BED} -b ${!t} -bedpe -type ospan -f 0.75 -s > ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bam}${t}.bedpe" >> $parafly_file
done
[ ! -f ${parafly_file}.completed ] && \
	ParaFly -c $parafly_file -CPU $CPU


# @ running bedpe S/AS counting
parafly_file=${INTERSECTOUTDIR}/countPE.para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	echo "awk 'BEGIN{OFS=\"\\t\"} {if (ARGIND==1) {ct[\$4]=0; ct_pos[\$4]=0; ct_neg[\$4]=0;} if (ARGIND==2) {ct[\$14]++; if (\$9==\$16){ct_pos[\$14]++} else {ct_neg[\$14]++}}}END{for (i in ct){ print i,ct[i],ct_pos[i],ct_neg[i] }}' ${!t} ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bam}${t}.bedpe  > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bam}${t}.bedpe.ct " >> $parafly_file && \
	echo "awk 'BEGIN{OFS=\"\\t\"} {if (ARGIND==1) {ct[\$4]=0; ct_pos[\$4]=0; ct_neg[\$4]=0;} if (ARGIND==2) {ct[\$14]++; if (\$9==\$16){ct_pos[\$14]++} else {ct_neg[\$14]++}}}END{for (i in ct){ print i,ct[i],ct_pos[i],ct_neg[i] }}' ${!t} ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bam}${t}.bedpe  > ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bam}${t}.bedpe.ct " >> $parafly_file
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

# @ running weblogo drawing
parafly_file=${LOGODIR}/weblogo.para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	echo "${RUN_PIPELINE_ADD}/bedPE2weblogo.sh ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bam}${t}.bedpe $GENOME_FA $5 ${LOGODIR}" >> $parafly_file && \
	echo "${RUN_PIPELINE_ADD}/bedPE2weblogo.sh ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bam}${t}.bedpe $GENOME_FA $5 ${LOGODIR}" >> $parafly_file
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

###########
#join pdfs#
###########
PDF_NAMES=""
for t in ${TARGETS[@]}
do \
	[ -f ${LOGODIR}/${UNIQUE_BED_NAME%bam}${t}.bedpe.reads.5end_60.information_content.pdf ]  && PDF_NAMES=${PDF_NAMES}" "${LOGODIR}/${UNIQUE_BED_NAME%bam}${t}.bedpe.reads.5end_60.information_content.pdf
	[ -f ${LOGODIR}/${MULTIP_BED_NAME%bam}${t}.bedpe.reads.5end_60.information_content.pdf ]  && PDF_NAMES=${PDF_NAMES}" "${LOGODIR}/${MULTIP_BED_NAME%bam}${t}.bedpe.reads.5end_60.information_content.pdf
done
OUTDIR=${INTERSECTOUTDIR%/*}
/home/wangw1/src/ghostscript-9.09/bin/gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${OUTDIR}/${PREFIX}.PE.pdf ${PDF_NAMES}

