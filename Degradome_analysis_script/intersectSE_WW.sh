#! /bin/bash -x

#####################
#variable assignment#
#####################
UNIQUE_BED=${1}
UNIQUE_BED_NAME=${1##*/}
MULTIP_BED=$2
MULTIP_BED_NAME=${2##*/}
OUT=$3
GENOME_FA=$5
[ ! -z $6 ] && CPU=$6 || CPU=8
[ ! -z $7 ] && PREFIX=$7 || PREFIX=${UNIQUE_BED_NAME%.uniq*}
INTERSECTOUTDIR=$8

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


# calculate NTM, concatenate UNIQ and MULT
[ ! -f ${UNIQUE_BED}.ntm ] && \
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,"1",$6}' $UNIQUE_BED > ${UNIQUE_BED}.ntm
	
[ ! -f ${UNIQUE_BED%unique.bed}all.bed.ntm ] && \
	bedtools groupby -i $MULTIP_BED -g 1,2,3,4,5,6 -c 4 -o distinct | \
	cut -f1-6 | \
	sort -k4,4 -k5,5 | \
	bedtools groupby -i - -g 4,5 -c 1,2,3,6 -o collapse,collapse,collapse,collapse | \
	awk 'BEGIN{OFS="\t"}{a=split($3,ar,","); split($4,br,",");  split($5,cr,","); split($6,dr,","); for(i=1;i<=a;i++){print ar[i],br[i],cr[i],$1,1.0/a,dr[i]}}' > ${MULTIP_BED}.ntm
	
[ ! -f ${UNIQUE_BED%unique.bed}all.bed.ntm ] && \
	cat ${UNIQUE_BED}.ntm ${MULTIP_BED}.ntm > ${UNIQUE_BED%unique.bed}all.bed.ntm
	
#collapse 5'end and sum the normalized read count up
[ ! -f ${UNIQUE_BED%unique.bed}all.bed.ntm.collapse ] && \
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,$4,$5}' ${UNIQUE_BED%unique.bed}all.bed.ntm | \
	bedtools sort -i - | \
	bedtools groupby -i - -g 1,2,3,4 -c 5,6 -o collapse,sum | \
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,$6,$4}' > ${UNIQUE_BED%unique.bed}all.bed.ntm.collapse
	
[ ! -f ${UNIQUE_BED}.ntm.collapse ] && \
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,$4,$5}' ${UNIQUE_BED}.ntm | \
	bedtools sort -i - | \
	bedtools groupby -i - -g 1,2,3,4 -c 5,6 -o collapse,sum | \
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,$6,$4}' > ${UNIQUE_BED}.ntm.collapse

#TAATAAGGCAATGATTGTTCAGTGCTTG	1	chrU:8868506-8868533(+)	antisense	repeats	suste.206	1	3
parafly_file=${INTERSECTOUTDIR}/intersect1.${RANDOM}.para && \
rm -rf $parafly_file
# processing data
for t in ${TARGETS[@]}
do \
	echo -ne "bedtools intersect -a ${UNIQUE_BED}.ntm.collapse -b ${!t} -f 0.99 -wb |bedtools sort -i - |bedtools groupby -i - -g 1,2,3,4,5,6 -c 7,8,9,10,11,12 -o collapse,collapse,collapse,collapse,collapse,collapse | \
		 awk 'BEGIN{OFS=\"\\\t\"}{a=split(\$7,ar,\",\");b=split(\$8,br,\",\"); c=split(\$9,cr,\",\");d=split(\$10,dr,\",\");e=split(\$11,er,\",\");f=split(\$12,fr,\",\"); for(i=1;i<=a;i++){print \$1,\$2,\$3,\$4,\$5,\$6,ar[i],br[i],cr[i],dr[i],er[i],fr[i],a} }'  \
		   > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME}.ntm.collapse.${t}.nta.bed && " >> $parafly_file ; 
	echo -ne "awk 'BEGIN{OFS=\"\\\t\"}{if(\$6==\$12){print \$4,\$5,\$1\":\"\$2\"-\"\$3\"(\"\$6\")\",\"sense\",\$10,\$10,\"1\",\$13}else{print \$4,\$5,\$1\":\"\$2\"-\"\$3\"(\"\$6\")\",\"antisense\",\$10,\$10,\"1\",\$13 } }' ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME}.ntm.collapse.${t}.nta.bed > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME}.ntm.collapse.${t}.nta.mapper2 &&" >> $parafly_file ;
	echo -e "gzip ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME}.ntm.collapse.${t}.nta.mapper2" >> $parafly_file ;
	
	echo -ne "bedtools intersect -a ${UNIQUE_BED%unique.bed}all.bed.ntm.collapse -b ${!t} -f 0.99 -wb |bedtools sort -i - |bedtools groupby -i - -g 1,2,3,4,5,6 -c 7,8,9,10,11,12 -o collapse,collapse,collapse,collapse,collapse,collapse | \
	     awk 'BEGIN{OFS=\"\\\t\"}{a=split(\$7,ar,\",\");b=split(\$8,br,\",\"); c=split(\$9,cr,\",\");d=split(\$10,dr,\",\");e=split(\$11,er,\",\");f=split(\$12,fr,\",\"); for(i=1;i<=a;i++){print \$1,\$2,\$3,\$4,\$5,\$6,ar[i],br[i],cr[i],dr[i],er[i],fr[i],a} }'\
		   > ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%unique.bed}all.bed.ntm.collapse.${t}.nta.bed && " >> $parafly_file ; 
	echo -ne "awk 'BEGIN{OFS=\"\\\t\"}{if(\$6==\$12){print \$4,\$5,\$1\":\"\$2\"-\"\$3\"(\"\$6\")\",\"sense\",\$10,\$10,\"1\",\$13}else{print \$4,\$5,\$1\":\"\$2\"-\"\$3\"(\"\$6\")\",\"antisense\",\$10,\$10,\"1\",\$13 } }' ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%unique.bed}all.bed.ntm.collapse.${t}.nta.bed >${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%unique.bed}all.bed.ntm.collapse.${t}.nta.mapper2 && " >> $parafly_file ;
	echo -e "gzip ${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%unique.bed}all.bed.ntm.collapse.${t}.nta.mapper2 " >> $parafly_file ;
done
# running ParaFly if no jobs has been ran (no .completed file) or it has ran but has some failed (has .failed_commands)
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU $CPU -failed_cmds $parafly_file.failed_commands
fi


