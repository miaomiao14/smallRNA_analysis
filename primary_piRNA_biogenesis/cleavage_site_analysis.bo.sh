CHROM=$HOME/piper/common/dm3/dm3.ChromInfo.txt
export PIPELINE_DIRECTORY=$HOME/piper
TOPN=${5:-10000}
BOTN=20000
CPU=8
[ ! -f ${1%bed2}piRNA.bed2 ] && awk '$3-$2>22 && $3-$2<30' $1 > ${1%bed2}piRNA.bed2
[ ! -f ${2%bed2}piRNA.bed2 ] && awk '$3-$2>22 && $3-$2<30' $2 > ${2%bed2}piRNA.bed2
##bed2Detail: this program calculate the ping-pong score for each genomic posisitons.
[ ! -f ${1%.ovary*}__${2%.ovary*}.bed2Detail ] && \
bed2Detail -a ${1%bed2}piRNA.bed2 -b ${2%bed2}piRNA.bed2 -c $CHROM | \
	awk 'BEGIN{OFS="\t"}{if ($6=="+"){print $1,$2,$2+1,$4/$5,1,$6,$8} else {print $1,$3-1,$3,$4/$5,1,$6,$8}}' | \
	gsort -k1,1 -k2,2 -k3,3 -k6,6 -T $PWD --parallel=$CPU | \
	bedtools groupby -i stdin -g 1,2,3,6 -c 4,7 -o sum,max | \
	gsort -k 6,6rn -T $PWD --parallel=$CPU > \
	${1%.ovary*}__${2%.ovary*}.bed2Detail
[ ! -f ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN} ] && awk 'BEGIN{OFS="\t"}{if ($4=="+") {print $1,$2,$3,$5,$6,$4} }' ${1%.ovary*}__${2%.ovary*}.bed2Detail | head -${TOPN} > ${1%.ovary*}__${2%.ovary*}.bed2Detail.Watson.top${TOPN}
[ ! -f ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN} ]  && awk 'BEGIN{OFS="\t"}{if ($4=="-") {print $1,$2,$3,$5,$6,$4} }' ${1%.ovary*}__${2%.ovary*}.bed2Detail | head -${TOPN} > ${1%.ovary*}__${2%.ovary*}.bed2Detail.Crick.top${TOPN}

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

