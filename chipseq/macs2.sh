#!/bin/sh


##assuming the input bam files are sorted
CHIPBAM=$1
InputBAM=$2

OUTDIR=$3
NAME=$4


macs2 callpeak -t $CHIPBAM -c $InputBAM -g dm -n $NAME -B --SPMR --outdir ${OUTDIR}
macs2 bdgcmp -t $NAME\_treat_pileup.bdg -c $NAME\_control_lambda.bdg -o ${OUTDIR}/$NAME.foldenrichment.bedGraph -m FE
macs2 bdgcmp -t $NAME\_treat_pileup.bdg -c $NAME\_control_lambda.bdg -o ${OUTDIR}/$NAME.logLR.bedGraph -m logLR -p 0.00001


if [ 0 -eq 0 ]
then
	d=`grep "^# d = " ${OUTDIR}/$NAME\_peaks.xls | cut -f4 -d " "` 
else
	d=0
fi

norm_treat=`grep "# total fragments in treatment:" ${OUTDIR}/$NAME\_peaks.xls | cut -f6 -d " "`
echo $norm_treat
norm_ctrl=`grep "# total fragments in control:" ${OUTDIR}/$NAME\_peaks.xls | cut -f6 -d " "`
echo $norm_ctrl


### bigWigs for browser###
bedtools bamtobed -i $CHIPBAM | awk -F "\t" -v d=$d '{ OFS="\t"; if($6=="-") { print $1,$3-d,$3 } else { print $1,$2,$2+d } }' | bedClip stdin /home/ww74w/pipeline/common/dm3.ChromInfo.txt tmp
sort -k1,1 -k2,2n tmp | bedItemOverlapCount null stdin -chromSize=/home/ww74w/pipeline/common/dm3.ChromInfo.txt | awk -F "\t" -v norm=$norm_treat '{OFS="\t"; print $1,$2,$3,$4*1000000/norm }' > tmp.bedGraph
bedGraphToBigWig tmp.bedGraph /home/ww74w/pipeline/common/dm3.ChromInfo.txt ${OUTDIR}/$NAME.treatment.bw

bedtools bamtobed -i $InputBAM | awk -F "\t" -v d=$d '{ OFS="\t"; if($6=="-") { print $1,$3-d,$3 } else { print $1,$2,$2+d } }' | bedClip stdin /home/ww74w/pipeline/common/dm3.ChromInfo.txt tmp
sort -k1,1 -k2,2n tmp | bedItemOverlapCount null stdin -chromSize=/home/ww74w/pipeline/common/dm3.ChromInfo.txt | awk -F "\t" -v norm=$norm_ctrl '{OFS="\t"; print $1,$2,$3,$4*1000000/norm }' > tmp.bedGraph
bedGraphToBigWig tmp.bedGraph /home/ww74w/pipeline/common/dm3.ChromInfo.txt ${OUTDIR}/$NAME.control.bw
