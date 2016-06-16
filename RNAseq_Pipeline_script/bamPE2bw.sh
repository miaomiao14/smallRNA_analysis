#!/usr/bin/env bash

export SPUCHROM=/diag/home/netashawang/cloud/spu_star_new/chrNameLength.txt

BAM=$1
STAROUT=$2
OUTDIR=$3 #output OUTDIR

#parafile=${OUTDIR}/${RANDOM}${RANDOM}.para
#touch $parafile

#for inbam in $BAM
#do

##Extract the number of unique mappers from Log.final.out of STAR
TC=`awk -F "\t" '/Uniquely mapped reads number/ {print $2}' $STAROUT`
echo $TC
let "NF=1000000/$TC"
echo $NF

[ -z "${prefix}" ] && prefix=${BAM%.bam}
##extract unique mappers
paraFile=${RANDOM}.para && \
samtools  view -b -q 255 $BAM | tee ${prefix}.unique.bam | bedtools bamtobed 
##bam -> bed
[ ! -f ${OUTDIR}/tmp.bed ] && \
bedtools bamtobed -i $inbam -bedpe -split | awk -F "\t" '{ OFS="\t"; if($6=="-") { print $1,$2,$3,".",0,"-" } else { print $1,$2,$3,".",0,"+" } }' | bedClip stdin $SPUCHROM ${OUTDIR}/tmp.bed && \
## pileup: bed -> bedGraph
#[ ! -f ${OUTDIR}/$prefix.bedGraph ] && \
#cut -f1-3 ${OUTDIR}/tmp.bed | sort -k1,1 -k2,2n | bedItemOverlapCount $gdb stdin -chromSize=$SPUCHROM | awk -F "\t" -v normFac=$normFac '{ OFS="\t"; print $1,$2,$3,$4*1000000/normFac; }' >${OUTDIR}/$prefix.bedGraph && \

## bedGraph -> bw
#[ ! -f ${OUTDIR}/$prefix.bw ] && \
#bedGraphToBigWig ${OUTDIR}/$prefix.bedGraph $SPUCHROM ${OUTDIR}/$prefix.bw

#done


