#! /bin/bash -x



DEG_PIPELINE_ADD=/home/ww74w/bin/
BEDPE=$1
BEDPE_NAME=${1##*/}
Genome=$2
STR=$3 # if STR==1, dUTR method, then take col4,5,6; otherwise take col1,2,3
SMRNA_PIPELINE_DIRECTORY=/home/ww74w/git/smallRNA_analysis/

##separate the piRNAs assigned to Aub, Ago3, or PIWI (from each IP) by length
##chrX	3302399	3302424	1	1	+	TTTTTTTTTTTTTTTTTTTTTAGGA

for i in {23..29}
do 
	awk -v l=$i 'BEGIN{OFS="\t"} len=$3-$2+1;if($len==$l) {print $0}' ${BEDPE_NAME}.$i
done


# doing by reads; 
# we don't use slopBed here because we do NOT want to resize according to the chromosome; since we want to make sure that all the sequences have the same length
# LOGODIR=$4
# [ -s $BEDPE ] && \
# if [[ $STR == 1 ]]
# then
# 	awk 'BEGIN{OFS="\t"}{print $4,$5,$6,1,1,$10}' $BEDPE > ${LOGODIR}/${BEDPE_NAME%bedpe}bed1
# else
# 	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,1,1,$9}' $BEDPE > ${LOGODIR}/${BEDPE_NAME%bedpe}bed1
# fi && \
# awk 'BEGIN{OFS="\t"} { if ($2>=30) { for (i=1;i<=$4;++i) { if ($6=="+") { print $1,$2-30,$2+31,$4,$5,$6 } else { print $1,$3-31,$3+30,$4,$5,$6 }}}}' ${LOGODIR}/${BEDPE_NAME%bedpe}bed1 > ${LOGODIR}/${BEDPE_NAME}.reads.5end_60 && \
# bedtools getfasta -fi $Genome -bed ${LOGODIR}/${BEDPE_NAME}.reads.5end_60 -fo stdout -s -name | tr 'tT' 'uU' > ${LOGODIR}/${BEDPE_NAME}.reads.5end_60.fa   && \
# $SMRNA_PIPELINE_DIRECTORY/bin/information_content.py  ${LOGODIR}/${BEDPE_NAME}.reads.5end_60.fa > ${LOGODIR}/${BEDPE_NAME}.reads.5end_60.information_content && \
# $SMRNA_PIPELINE_DIRECTORY/bin/Rscript $SMRNA_PIPELINE_DIRECTORY/bin/draw_percentage.R ${LOGODIR}/${BEDPE_NAME}.reads.5end_60.information_content ${LOGODIR}/${BEDPE_NAME}.reads.5end_60.information_content && \
# rm -rf ${LOGODIR}/${BEDPE_NAME%bedpe}bed1   ${LOGODIR}/${BEDPE_NAME}.reads.5end_60   ${LOGODIR}/${BEDPE_NAME}.reads.5end_60.fa