#!/bin/bash 

CHRSIZE=$1
OUTDIR=$2
prefix=$3
LOG=$4

#calculate the genome coverage of each annotation feature, print it to the lOG file
rm $LOG
awk '{OFS="\t"}{print $1,0,$2-1}' $CHRSIZE >${OUTDIR}/${prefix}.chrsize.bed
echo "FirstExonCoverage GenomeSize FirstExonCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.firstexon.bed -b ${OUTDIR}/${prefix}.chrsize.bed > ${OUTDIR}/${prefix}.sorted.clipped.firstexon.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.firstexon.bed.coverage >> $LOG
echo "LastExonCoverage GenomeSize LastExonCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.lastexon.bed -b ${OUTDIR}/${prefix}.chrsize.bed >${OUTDIR}/${prefix}.sorted.clipped.lastexon.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.lastexon.bed.coverage  >> $LOG
echo "MiddleExonCoverage GenomeSize MiddleExonCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.middleexon.bed -b ${OUTDIR}/${prefix}.chrsize.bed >${OUTDIR}/${prefix}.sorted.clipped.middleexon.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.middleexon.bed.coverage  >> $LOG
echo "SingleExonCoverage GenomeSize SingleExonCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.singleexon.bed -b ${OUTDIR}/${prefix}.chrsize.bed >${OUTDIR}/${prefix}.sorted.clipped.singleexon.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.singleexon.bed.coverage >> $LOG
echo "AllExonCoverage GenomeSize AllExonCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.allexon.bed -b ${OUTDIR}/${prefix}.chrsize.bed >${OUTDIR}/${prefix}.sorted.clipped.allexon.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.allexon.bed.coverage >> $LOG
echo "AllIntronCoverage GenomeSize AllIntronCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.intron.bed -b ${OUTDIR}/${prefix}.chrsize.bed >${OUTDIR}/${prefix}.sorted.clipped.intron.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.intron.bed.coverage >> $LOG
echo "ExonIntronJunCoverage GenomeSize ExonIntronJunCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.ExonIntronJun.slopped.bed -b ${OUTDIR}/${prefix}.chrsize.bed >${OUTDIR}/${prefix}.sorted.clipped.ExonIntronJun.slopped.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.ExonIntronJun.slopped.bed.coverage >> $LOG
echo "IntronExonJunCoverage GenomeSize IntronExonJunCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.IntronExonJun.slopped.bed -b ${OUTDIR}/${prefix}.chrsize.bed >${OUTDIR}/${prefix}.sorted.clipped.IntronExonJun.slopped.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.IntronExonJun.slopped.bed.coverage >> $LOG
echo "ExonExonJunCoverage GenomeSize ExonExonJunCoverageFraction" >> $LOG
bedtools coverage -a ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.bed -b ${OUTDIR}/${prefix}.chrsize.bed >${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.bed.coverage
awk '{a=a+$5;b=b+$6;}END{print a,b,a/b}' ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.bed.coverage >> $LOG


