CHROM=dm3.ChromInfo.txt
TOPN=${5:-10000}
BOTN=20000
CPU=8
echo "taking the top $TOPN lines of input file and convert to bed"
awk 'BEGIN{OFS="\t"}{if ($4=="+") {printf "%s\t%d\t%d\t%.2f\t%d\t%s\n", $1,$2,$3,$5,$6,$4} }' $1 | head -${TOPN} > ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}
awk 'BEGIN{OFS="\t"}{if ($4=="-") {printf "%s\t%d\t%d\t%.2f\t%d\t%s\n", $1,$2,$3,$5,$6,$4} }' $1 | head -${TOPN} >  ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}
echo "use bwtool to aggregate the sigal"
bwtool agg 100:100 -starts ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN} \
	 $3 \
	 ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.txt
bwtool agg 100:100 -starts ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN} \
	 $4 \
	 ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.txt
paste ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.txt  ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.txt | awk 'BEGIN{OFS="\t"}{if ($1>0) --$1; print $1,$2-$4}' > ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}.txt && rm -rf ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.txt  ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.txt
Rscript --slave piPipes_draw_aggregate.R ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}.txt ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*} ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}

bwtool matrix 100:100 -starts ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN} \
	 $3 \
	 ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.matrix1 && \
sed -e 's/NA/0/g' ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.matrix1 > ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.matrix
bwtool matrix 100:100 -starts ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN} \
	 $4 \
	 ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.matrix1 && \
sed -e 's/NA/0/g' -e 's/-//g' ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.matrix1 > ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.matrix
cat ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.matrix ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.matrix > ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}.matrix && \
	awk 'BEGIN{OFS="\t"}{s=0; for (i=1;i<=NF;++i) s+=$i; print $0"\t"s}' ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}.matrix | sort -k201,201rn | cut -f1-200 > ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}.matrix1 && \
	mv ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}.matrix1 ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}.matrix && \
	Rscript --slave drawHP.R ${1%.ovary*}__${2%.ovary*}.bed2Detail.top${TOPN}_agg_${3%.ovary*}.matrix && \
	rm -rf ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.matrix1 ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}_agg_${3%.ovary*}.matrix ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.matrix1 ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}_agg_${3%.ovary*}.matrix

# [ ! -f ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.bot${BOTN} ] && awk 'BEGIN{OFS="\t"}{if ($4=="+") {print $1,$2,$3,$5,$6,$4} }' ${1%.ovary*}__${2%.ovary*}.bed2Detail | tail -${BOTN} > ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.bot${BOTN}
# [ ! -f ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.bot${BOTN} ]  && awk 'BEGIN{OFS="\t"}{if ($4=="-") {print $1,$2,$3,$5,$6,$4} }' ${1%.ovary*}__${2%.ovary*}.bed2Detail | tail -${BOTN} > ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.bot${BOTN}
#
# bwtool agg 100:100 -starts ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.bot${BOTN} \
# 	 $3 \
# 	 ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.bot${BOTN}_agg_${3%.ovary*}.txt
#
# bwtool agg 100:100 -starts ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.bot${BOTN} \
# 	 $4 \
# 	 ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.bot${BOTN}_agg_${3%.ovary*}.txt
# paste ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.bot${BOTN}_agg_${3%.ovary*}.txt  ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.bot${BOTN}_agg_${3%.ovary*}.txt | awk 'BEGIN{OFS="\t"}{if ($1>0) --$1; print $1,$2-$4}' > ${1%.ovary*}__${2%.ovary*}.bed2Detail.bot${BOTN}_agg_${3%.ovary*}.txt && rm -rf ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.bot${BOTN}_agg_${3%.ovary*}.txt  ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.bot${BOTN}_agg_${3%.ovary*}.txt
# Rscript --slave piPipes_draw_aggregate.R ${1%.ovary*}__${2%.ovary*}.bed2Detail.bot${BOTN}_agg_${3%.ovary*}.txt ${1%.ovary*}__${2%.ovary*}.bed2Detail.bot${BOTN}_agg_${3%.ovary*} ${1%.ovary*}__${2%.ovary*}.bed2Detail.bot${BOTN}_agg_${3%.ovary*}
