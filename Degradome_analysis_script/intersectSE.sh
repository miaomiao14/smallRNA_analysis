#! /bin/bash -x

#####################
#variable assignment#
#####################
RUN_PIPELINE_ADD=/home/wangw1/git/PE_RNA_Pipeline
UNIQUE_BED=${1}
UNIQUE_BED_NAME=${1##*/}
MULTIP_BED=$2
MULTIP_BED_NAME=${2##*/}
OUT=$3
GENOME_FA=$5
[ ! -z $6 ] && CPU=$6 || CPU=8
[ ! -z $7 ] && PREFIX=$7 || PREFIX=${UNIQUE_BED_NAME%.uniq*}
INTERSECTOUTDIR=$8
LOGODIR=$9
#test
echo $INTERSECTOUTDIR
echo $LOGODIR
case ${4} in
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
;;
*)
	echo "unknown orgnanism... currently only mouse/fly is supported"
	exit 2
;;
esac

##############
#print header#
##############
echo -ne "Sample\tTotal_Unique_Reads\t" > $OUT;
for t in ${TARGETS[@]}
do \
	echo -ne "${t}_uniq_reads\t${t}_sense_uniq_reads\t${t}_antisense_uniq_reads\t" >> $OUT;
done
echo -ne "\n" >> $OUT;

##################
#calculate counts#
##################
echo -ne "${UNIQUE_BED%%.bed*}\t" >> $OUT;
TOTAL_UNIQ_READS=`wc -l $UNIQUE_BED | cut -f1 -d' '`
echo -ne "$TOTAL_UNIQ_READS\t" >> $OUT;

parafly_file=${INTERSECTOUTDIR}/intersect1.para && \
rm -rf $parafly_file
# processing data
for t in ${TARGETS[@]}
do \
	echo "bedtools intersect -wa -u -a $UNIQUE_BED -b  ${!t} > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.bed" >> $parafly_file ; 
	echo "bedtools intersect -wa -u -a $MULTIP_BED -b  ${!t} > ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.bed" >> $parafly_file ; 
done
# running ParaFly if no jobs has been ran (no .completed file) or it has ran but has some failed (has .failed_commands)
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

parafly_file=${INTERSECTOUTDIR}/intersect2.para && \
rm -rf $parafly_file
# processing data
for t in ${TARGETS[@]}
do \
	echo "[ -s ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.bed ] && bedtools intersect -wa -u -a ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.bed -b  ${!t} -s > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.S.bed  " >> $parafly_file ;
	echo "[ -s ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.bed ] && bedtools intersect -wa -u -a ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.bed -b  ${!t} -S > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.AS.bed " >> $parafly_file ;
	echo "[ -s ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.bed ] && bedtools intersect -wa -u -a ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.bed -b  ${!t} -s > ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.S.bed  " >> $parafly_file ;
	echo "[ -s ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.bed ] && bedtools intersect -wa -u -a ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.bed -b  ${!t} -S > ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.AS.bed " >> $parafly_file ;
done
# running ParaFly if no jobs has been ran (no .completed file) or it has ran but has some failed (has .failed_commands)
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi

##################
#calculate counts#
##################
for t in ${TARGETS[@]}
do \
	TOTAL_READS=`wc -l ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.bed | cut -f1 -d' '`
	SENSE_READS=`wc -l ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.S.bed | cut -f1 -d' '`
	ANTISENSE_READS=`wc -l ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.AS.bed | cut -f1 -d' '`
	echo -ne $TOTAL_READS"\t"$SENSE_READS"\t"$ANTISENSE_READS"\t" >> $OUT;
done

# transpose the table
ruby -lane 'BEGIN{$arr=[]}; $arr.concat([$F]); END{$arr.transpose.each{ |a| puts a.join ("\t") } }' -F"\t" $OUT > ${OUT}.t  && \
mv ${OUT}.t ${OUT} 

# @ running weblogo drawing
parafly_file=${LOGODIR}/weblogo.para && \
rm -rf $parafly_file
for t in ${TARGETS[@]}
do \
	echo "${RUN_PIPELINE_ADD}/bedSE2weblogo.sh ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.S.bed  $GENOME_FA ${LOGODIR}" >> $parafly_file && \
	echo "${RUN_PIPELINE_ADD}/bedSE2weblogo.sh ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bed}${t}.AS.bed $GENOME_FA ${LOGODIR}" >> $parafly_file && \
	echo "${RUN_PIPELINE_ADD}/bedSE2weblogo.sh ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.S.bed  $GENOME_FA ${LOGODIR}" >> $parafly_file && \
	echo "${RUN_PIPELINE_ADD}/bedSE2weblogo.sh ${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bed}${t}.AS.bed $GENOME_FA ${LOGODIR}" >> $parafly_file
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
	[ -f ${LOGODIR}/${UNIQUE_BED_NAME%bed}${t}.S.bed.reads.5end_60.information_content.pdf ]  && PDF_NAMES=${PDF_NAMES}" "${LOGODIR}/${UNIQUE_BED_NAME%bed}${t}.S.bed.reads.5end_60.information_content.pdf
	[ -f ${LOGODIR}/${UNIQUE_BED_NAME%bed}${t}.AS.bed.reads.5end_60.information_content.pdf ] && PDF_NAMES=${PDF_NAMES}" "${LOGODIR}/${UNIQUE_BED_NAME%bed}${t}.AS.bed.reads.5end_60.information_content.pdf
	[ -f ${LOGODIR}/${MULTIP_BED_NAME%bed}${t}.S.bed.reads.5end_60.information_content.pdf ]  && PDF_NAMES=${PDF_NAMES}" "${LOGODIR}/${MULTIP_BED_NAME%bed}${t}.S.bed.reads.5end_60.information_content.pdf
	[ -f ${LOGODIR}/${MULTIP_BED_NAME%bed}${t}.AS.bed.reads.5end_60.information_content.pdf ] && PDF_NAMES=${PDF_NAMES}" "${LOGODIR}/${MULTIP_BED_NAME%bed}${t}.AS.bed.reads.5end_60.information_content.pdf
done
OUTDIR=${INTERSECTOUTDIR%/*}
/home/wangw1/src/ghostscript-9.09/bin/gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${OUTDIR}/${PREFIX}.SE.pdf ${PDF_NAMES}



