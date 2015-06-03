#! /usr/bin/env bash

BED=$1
BED_NAME=${1##*/}
Genome=$2
LOGODIR=$3
SMRNA_PIPELINE_DIRECTORY=/home/hanb/nearline/small_RNA_Pipeline

# doing by reads; we don't use slopBed here because we do NOT want to resize according to the chromosome; since we want to make sure that all the sequences have the same length
[ -s $BED ] && \
	awk 'BEGIN{OFS="\t"} { if ($6=="+") { print $1,$2-30,$2+31,$4,$5,$6 } else { print $1,$3-31,$3+30,$4,$5,$6 }}' $BED | \
	awk '$2>0 && $3>0' > ${LOGODIR}/${BED_NAME}.reads.5end_60 && \
	bedtools getfasta -fi $Genome -bed ${LOGODIR}/${BED_NAME}.reads.5end_60 -fo stdout -s -name | tr 'tT' 'uU' > \
		${LOGODIR}/${BED_NAME}.reads.5end_60.fa  && \
	$SMRNA_PIPELINE_DIRECTORY/bin/information_content.py  ${LOGODIR}/${BED_NAME}.reads.5end_60.fa > ${LOGODIR}/${BED_NAME}.reads.5end_60.information_content && \
	$SMRNA_PIPELINE_DIRECTORY/bin/Rscript $SMRNA_PIPELINE_DIRECTORY/bin/draw_percentage.R ${LOGODIR}/${BED_NAME}.reads.5end_60.information_content ${LOGODIR}/${BED_NAME}.reads.5end_60.information_content
	rm -rf ${LOGODIR}/${BED_NAME}.reads.5end_60 ${LOGODIR}/${BED_NAME}.reads.5end_60.fa