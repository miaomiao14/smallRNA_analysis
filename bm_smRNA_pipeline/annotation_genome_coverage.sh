#!/bin/sh 

CHRSIZE=$1
INDIR=$2
OUTDIR=$3
[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
prefix="bm"
LOG=${OUTDIR}/coverage.log

#calculate the genome coverage of each annotation feature, print it to the lOG file
rm $LOG
awk '{OFS="\t"}{print $1,0,$2-1}' $CHRSIZE >${OUTDIR}/${prefix}.chrsize.bed

for i in `ls ${INDIR}/silkworm.*.bed`
do
bedtools coverage -a ${i} -b ${OUTDIR}/${prefix}.chrsize.bed > ${OUTDIR}/${i##*/}.coverage && \
echo -e "${i##*/}\tGenomeSize\tCoverageFraction" >> $LOG && \
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${i}.coverage >> $LOG

done