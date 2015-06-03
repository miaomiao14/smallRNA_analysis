#! /bin/bash -x
case ${3} in
## mouse specific intersectBed target files
"mouse")
	FOLDER=/home/hanb/nearline/small_RNA_Pipeline/common_files/mm9/UCSC_BEDS
	# piRNA clusters annoated in Li, et al, Mol Cell, 2013
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


#UNIQNORMBED=$1
#UNIQNORMBED_NAME=${1##*/}
MULTNORMBED=$1
MULTNORMBED_NAME=${1##*/}
#annotation_bed=$2
OUTDIR=$2
READ_LEN=$4
#normbed file are from bamtobed and the 5th column of the normbed file is 1/NTM

#bedtools bed12tobed6 -i $annotation_bed | sort -u | sort -k1,1 -k2,2n > annotation_bed 

# 1. split splicing alignment
# 2. intersect with the gene annotation
# 3. remove reduant records caused by isoform
# 4. sum the total overlap for the splicing alignment
# 5. check if the two ends are both fully covered in the same feature

parafly_file=${OUTDIR}/bed6intersecttobed6.${RANDOM}.para && \
for t in ${TARGETSBED6[@]}
do \
	#bed12tobed6 split the splicing alignment; each splicing part of one read has to be in the spliced annotation feature, then the sum of intersect 
	#bedtools bed12tobed6 -i $MULTNORMBED | bedtools intersect -a stdin -b annotation_bed -wo | bedtools groupby -g 1,2,3,4,5,6,10,12 -c 13 -o max | bedtools groupby -g 4,5,6,7,8 -c 9 -o sum | sed 's/\/[1|2]//;' | sort -k1,1 > tmp
	
	echo -ne "cat ${MULTNORMBED} |awk '{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$7\"/1\",\$8,\$9\"\\\n\"\$4,\$5,\$6,\$7\"/2\",\$8,\$10}' | bedtools intersect -a stdin -b ${!t} -wo | bedtools groupby -g 1,2,3,4,5,6,10,12 -c 13 -o max | bedtools groupby -g 4,5,6,7,8 -c 9 -o sum | sed 's/\/[1|2]//;' | sort -k1,1 > ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.tmp && " >> $parafly_file
	# sense; remove the reads annotated with multiple features(if a read overlap with multiple features, it's not included in the final analysis) and the alignment length has to be readLenth
	echo -ne "awk '{OFS=\"\\\t\"; if(\$3==\$5) print }' ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.tmp | bedtools groupby -g 1,2,3, -c 4,5,6 -o distinct,distinct,distinct | awk -v len=$READ_LEN '{OFS=\"\\\t\"; }{if(!/,/ && \$6==len) print \$4,\$2}' | sort -k1,1 > ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.sense && ">> $parafly_file  
	# sense, allmap
	echo -ne "bedtools groupby -i ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.sense -g 1 -c 2 -o sum > ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.sense.fragment.counts && ">> $parafly_file
	# sense, uniqmap
	echo -ne "awk '{OFS=\"\\\t\"; if(\$2==1) print }' ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.sense | bedtools groupby -g 1 -c 2 -o sum > ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.uniqmap.sense.fragment.counts && ">> $parafly_file

	# antisense
	echo -ne "awk '{OFS=\"\\\t\"; if(\$3!=\$5) print }' ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.tmp | bedtools groupby -g 1,2,3, -c 4,5,6 -o distinct,distinct,distinct | awk -v len=$READ_LEN '{OFS=\"\\\t\";}{ if(!/,/ && \$6==len) print \$4,\$2}' | sort -k1,1 > ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.antisense && ">> $parafly_file 
	# antisense, allmap
	echo -ne "bedtools groupby -i ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.antisense -g 1 -c 2 -o sum > ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.antisense.fragment.counts && ">> $parafly_file
	# antisense, uniqmap
	echo -e "awk '{OFS=\"\\\t\"; if(\$2==1) print }' ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.antisense | bedtools groupby -g 1 -c 2 -o sum > ${OUTDIR}/${MULTNORMBED_NAME%bedpe}${t}.uniqmap.antisense.fragment.counts ">> $parafly_file
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU 4 -failed_cmds $parafly_file.failed_commands
fi
#rm ${OUTDIR}/*.tmp
#rm ${OUTDIR}/*.sense
#rm ${OUTDIR}/*.antisense