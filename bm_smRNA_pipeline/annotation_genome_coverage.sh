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

for i in `ls ${INDIR}/*.bed`
do
	name=${i##*/}
bedtools coverage -a ${i} -b ${OUTDIR}/${prefix}.chrsize.bed > ${OUTDIR}/${name}.coverage && \
echo -e "${i##*/}\tGenomeSize\tCoverageFraction" >> $LOG && \
awk '{OFS="\t"}{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${i}.coverage >> $LOG

done