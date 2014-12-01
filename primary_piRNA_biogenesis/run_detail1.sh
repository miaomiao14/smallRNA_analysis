#para1: input bed file(anchor point): chr start end strand reads ppscore
#para2: _ 
#para3: bigWig files for piRNAs (watson )
#para4: bigWig files for piRNAs (crick)
#para5: can be ignored

pipeline_Dir=${HOME}/git/smallRNA_analysis/primary_piRNA_biogenesis

CHROM=${pipeline_Dir}/dm3.ChromInfo.txt
TOPN=${5:-10000}
OUTDIR=$6
BOTN=20000
CPU=8

outfileName1=${1##*/}
#outfileName2=${2##*/}
outfileName3=${3##*/}
#outfileName3=${outfileName3/UA_VA.ppseq.txt./}
outfileName3=${outfileName3%.ovary*}
#namePrefix=${1%.ovary*}__
namePrefix=${outfileName1/UA_VA.ppseq.txt./}__

echo "taking the top $TOPN lines of input file and convert to bed"
awk 'BEGIN{OFS="\t"}{if ($4=="+") {printf "%s\t%d\t%d\t%.2f\t%d\t%s\n", $1,$2,$3,$5,$6,$4} }' $1 | head -${TOPN} > ${namePrefix}.bed2Detail.Watson.top${TOPN}
awk 'BEGIN{OFS="\t"}{if ($4=="-") {printf "%s\t%d\t%d\t%.2f\t%d\t%s\n", $1,$2,$3,$5,$6,$4} }' $1 | head -${TOPN} >  ${namePrefix}.bed2Detail.Crick.top${TOPN}
echo "use bwtool to aggregate the sigal"
#waston strand


bwtool agg 100:100 -starts ${namePrefix}.bed2Detail.Watson.top${TOPN} \
	 $3 \
	 ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.txt
#crick strand
bwtool agg 100:100 -starts ${namePrefix}.bed2Detail.Crick.top${TOPN} \
	 $4 \
	 ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.txt
#merge waston and crick
paste ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.txt  ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.txt | awk 'BEGIN{OFS="\t"}{if ($1>0) --$1; print $1,$2-$4}' > ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}.txt && rm -rf ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.txt  ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.txt
#
Rscript --slave ${pipeline_Dir}/piPipes_draw_aggregate.R ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}.txt ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*} ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}

bwtool matrix 100:100 -starts ${namePrefix}.bed2Detail.Watson.top${TOPN} \
	 $3 \
	 ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.matrix1 && \
sed -e 's/NA/0/g' ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.matrix1 > ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.matrix
bwtool matrix 100:100 -starts ${namePrefix}.bed2Detail.Crick.top${TOPN} \
	 $4 \
	 ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.matrix1 && \
sed -e 's/NA/0/g' -e 's/-//g' ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.matrix1 > ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.matrix
cat ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.matrix ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.matrix > ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}.matrix && \
	awk 'BEGIN{OFS="\t"}{s=0; for (i=1;i<=NF;++i) s+=$i; print $0"\t"s}' ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}.matrix | sort -k201,201rn | cut -f1-200 > ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}.matrix1 && \
	mv ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}.matrix1 ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}.matrix && \
	Rscript --slave ${pipeline_Dir}/drawHP.R ${OUTDIR}/${namePrefix}.bed2Detail.top${TOPN}_agg_${outfileName3%.ovary*}.matrix && \
	rm -rf ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.matrix1 ${OUTDIR}/${namePrefix}.bed2Detail.Watson.top${TOPN}_agg_${outfileName3%.ovary*}.matrix ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.matrix1 ${OUTDIR}/${namePrefix}.bed2Detail.Crick.top${TOPN}_agg_${outfileName3%.ovary*}.matrix

# [ ! -f ${namePrefix}.bed2Detail.Watson.bot${BOTN} ] && awk 'BEGIN{OFS="\t"}{if ($4=="+") {print $1,$2,$3,$5,$6,$4} }' ${namePrefix}.bed2Detail | tail -${BOTN} > ${namePrefix}.bed2Detail.Watson.bot${BOTN}
# [ ! -f ${namePrefix}.bed2Detail.Crick.bot${BOTN} ]  && awk 'BEGIN{OFS="\t"}{if ($4=="-") {print $1,$2,$3,$5,$6,$4} }' ${namePrefix}.bed2Detail | tail -${BOTN} > ${namePrefix}.bed2Detail.Crick.bot${BOTN}
#
# bwtool agg 100:100 -starts ${namePrefix}.bed2Detail.Watson.bot${BOTN} \
# 	 $3 \
# 	 ${namePrefix}.bed2Detail.Watson.bot${BOTN}_agg_${outfileName3%.ovary*}.txt
#
# bwtool agg 100:100 -starts ${namePrefix}.bed2Detail.Crick.bot${BOTN} \
# 	 $4 \
# 	 ${namePrefix}.bed2Detail.Crick.bot${BOTN}_agg_${outfileName3%.ovary*}.txt
# paste ${namePrefix}.bed2Detail.Watson.bot${BOTN}_agg_${outfileName3%.ovary*}.txt  ${namePrefix}.bed2Detail.Crick.bot${BOTN}_agg_${outfileName3%.ovary*}.txt | awk 'BEGIN{OFS="\t"}{if ($1>0) --$1; print $1,$2-$4}' > ${namePrefix}.bed2Detail.bot${BOTN}_agg_${outfileName3%.ovary*}.txt && rm -rf ${namePrefix}.bed2Detail.Watson.bot${BOTN}_agg_${outfileName3%.ovary*}.txt  ${namePrefix}.bed2Detail.Crick.bot${BOTN}_agg_${outfileName3%.ovary*}.txt
# Rscript --slave piPipes_draw_aggregate.R ${namePrefix}.bed2Detail.bot${BOTN}_agg_${outfileName3%.ovary*}.txt ${namePrefix}.bed2Detail.bot${BOTN}_agg_${outfileName3%.ovary*} ${namePrefix}.bed2Detail.bot${BOTN}_agg_${outfileName3%.ovary*}
