#!/bin/sh

INPUT=$1
OUTDIR=$2
CUTOFF=$3

##filter peaks with q value >4

filename=${INPUT##*/}
awk -v c=$CUTOFF 'BEGIN{OFS="\t"}{if($0 !~ /#/){ if($7 >=2 && $8 >= c) print $0 } }' $INPUT > ${OUTDIR}/$filename.qvalue.$CUTOFF.bed

##2231 peaks left in wildtype
##overlap with euchromatin and heterochromatin

hetschr=/home/ww74w/pipeline/common/dm3.hetChromInfo
#1231 peaks
bedtools intersect -a ${OUTDIR}/$filename.qvalue.$CUTOFF.bed -b $hetschr -wa -f 0.5 >${OUTDIR}/$filename.qvalue.$CUTOFF.heteroChromatin.bed

#1000 peaks
bedtools intersect -a ${OUTDIR}/$filename.qvalue.$CUTOFF.bed -b $hetschr -wa -f 0.5 -v >${OUTDIR}/$filename.qvalue.$CUTOFF.EuChromatin.bed