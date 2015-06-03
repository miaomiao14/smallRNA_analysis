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
	FOLDER=/home/hanb/nearline/small_RNA_Pipeline/common_files/dm3/UCSC_BEDS
	
	FLY_PIRNA_CLUSTER=$FOLDER/Brennecke.pirnaCluster.gtf 
	FLY_flyBase_GENE=$FOLDER/UCSC.flyBase.Genes.gtf
	FLY_flyBase_EXON=$FOLDER/UCSC.flyBase.Exon.gtf
	FLY_flyBase_INTRON=$FOLDER/UCSC.flyBase.Intron.gtf
	FLY_flyBase_5UTR=$FOLDER/UCSC.flyBase.5UTR.gtf
	FLY_flyBase_3UTR=$FOLDER/UCSC.flyBase.3UTR.gtf
	FLY_REPEATMASKER_DNA=$FOLDER/repeat_mask.gtf.DNA
	FLY_REPEATMASKER_LINE=$FOLDER/repeat_mask.gtf.LINE
	FLY_REPEATMASKER_LTR=$FOLDER/repeat_mask.gtf.LTR
	FLY_REPEATMASKER_RNA=$FOLDER/repeat_mask.gtf.RNA
	FLY_REPEATMASKER_rRNA=$FOLDER/repeat_mask.gtf.rRNA
	FLY_REPEATMASKER_Satellite=$FOLDER/repeat_mask.gtf.Satellite
	FLY_REPEATMASKER_Simple_repeat=$FOLDER/repeat_mask.gtf.Simple_repeat
	FLY_REPEATMASKER_Unknown=$FOLDER/repeat_mask.gtf.Unknown
	
	FLY_TRANSPOSON_ALL=$FOLDER/transposon.gtf2
	FLY_TRANSPOSON_GROUP1=$FOLDER/Zamore.group1.gtf
	FLY_TRANSPOSON_GROUP2=$FOLDER/Zamore.group2.gtf
	FLY_TRANSPOSON_GROUP3=$FOLDER/Zamore.group3.gtf
	FLY_TRANSPOSON_GROUP0=$FOLDER/Zamore.group0.gtf
	
	declare -a TARGETS=("FLY_PIRNA_CLUSTER" \
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

BAM=$1
BAM_NAME=${BAM##*/}
CAGELIKE=$2
OUTPUTDIR=$4
CPU=$5
case $CAGELIKE in
	1)
	LIB_TYPE="fr-firststrand"
	;;
	0)
	LIB_TYPE="fr-secondstrand"
	;;
esac
parafly_file=${OUTPUTDIR}/cuff.${RANDOM}.para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	OUTDIR=${OUTPUTDIR}/${t}
	[ ! -f ${OUTDIR} ] && mkdir -p ${OUTDIR}
	echo -ne "
cufflinks -v --no-update-check \
	--library-type ${LIB_TYPE} \
	-p ${CPU} \
	-g ${!t} \
	-G ${!t} \
	-o ${OUTDIR} \
	$BAM \
	-u \
	-j 0.2 \
	--min-frags-per-transfrag 40 \
	--overlap-radius 100 \
	2&> ${OUTDIR}/cufflinks.log && " >> $parafly_file
echo -ne " mv ${OUTDIR}/transcripts.gtf ${OUTPUTDIR}/${BAM_NAME%bam}${t}_cufflink_transcripts.gtf && " >> $parafly_file
echo -ne " mv ${OUTDIR}/genes.fpkm_tracking ${OUTPUTDIR}/${BAM_NAME%bam}${t}_cufflink_genes.fpkm_tracking && " >> $parafly_file
echo -e "mv ${OUTDIR}/isoforms.fpkm_tracking ${OUTPUTDIR}/${BAM_NAME%bam}${t}_cufflink_isoforms.fpkm_tracking " >> $parafly_file
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU 3 -failed_cmds $parafly_file.failed_commands
fi

