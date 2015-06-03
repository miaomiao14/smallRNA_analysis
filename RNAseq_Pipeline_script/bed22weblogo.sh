#! /usr/bin/env bash

BEDPE=$1
Genome=$2

# doing by reads; we don't use slopBed here because we do NOT want to resize according to the chromosome; since we want to make sure that all the sequences have the same length
[ -s $BEDPE ] && \
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,1,1,}' $BEDPE > ${BEDPE%bedpe}bed1 && \
awk 'BEGIN{OFS="\t"} { if ($2>=30) { for (i=1;i<=$4;++i) { if ($6=="+") { print $1,$2-30,$2+31,$4,$5,$6 } else { print $1,$3-31,$3+30,$4,$5,$6 }}}}' ${BEDPE%bedpe}bed1 > $BEDPE.reads.5end_60 && \
bedtools getfasta -fi $Genome -bed $BEDPE.reads.5end_60 -fo stdout -s -name | tr 'tT' 'uU' > $BEDPE.reads.5end_60.fa   && \
weblogo  --composition none -F PDF -A rna -W 15 -n 62 -c classic -S 1 < $BEDPE.reads.5end_60.fa > $BEDPE.reads.5end_60.seqlogo2.pdf && \
information_content.py  $BEDPE.reads.5end_60.fa > $BEDPE.reads.5end_60.information_content && \
Rscript $PIPELINE_DIRECTORY/bin/draw_information_content.R $BEDPE.reads.5end_60.information_content $BEDPE.reads.5end_60.information_content && \

#awk 'BEGIN{OFS="\t"} { if ($2>=30) { if ($6=="+") { print $1,$2-30,$2+31,$4,$5,$6 } else { print $1,$3-31,$3+30,$4,$5,$6 }}}' $BEDPE > $BEDPE.species.5end_60 && \
#bedtools getfasta -fi $Genome -bed $BEDPE.species.5end_60 -fo stdout -s -name | tr 'tT' 'uU' > $BEDPE.species.5end_60.fa && \
#weblogo  --composition none -F PDF -A rna -W 15 -n 62 -c classic -S 2 < $BEDPE.species.5end_60.fa > $BEDPE.species.5end_60.seqlogo.pdf && \
#weblogo  --composition none -F PDF -A rna -W 15 -n 62 -c classic -S 1 < $BEDPE.species.5end_60.fa > $BEDPE.species.5end_60.seqlogo2.pdf && \
#information_content.py  $BEDPE.species.5end_60.fa > $BEDPE.species.5end_60.information_content && \
#Rscript $PIPELINE_DIRECTORY/bin/draw_information_content.R $BEDPE.species.5end_60.information_content $BEDPE.species.5end_60.information_content && \
#awk 'BEGIN{OFS="\t"} { if ($2>=30) { if ($6=="+") { print $1,$2-30,$2+31,$4,$5,$6 } else { print $1,$3-31,$3+30,$4,$5,$6 }}}' $BEDPE | awk '{a[$0]=1}END{for (b in a) {print b}}' > $BEDPE.species.no3Hetero.5end_60 && \
#bedtools getfasta -fi $Genome -bed $BEDPE.species.no3Hetero.5end_60 -fo stdout -s -name | tr 'tT' 'uU' > $BEDPE.species.no3Hetero.5end_60.fa && \
#weblogo  --composition none -F PDF -A rna -W 15 -n 62 -c classic -S 2 < $BEDPE.species.no3Hetero.5end_60.fa > $BEDPE.species.no3Hetero.5end_60.seqlogo.pdf && \
#weblogo  --composition none -F PDF -A rna -W 15 -n 62 -c classic -S 1 < $BEDPE.species.no3Hetero.5end_60.fa > $BEDPE.species.no3Hetero.5end_60.seqlogo2.pdf && \
#information_content.py $BEDPE.species.no3Hetero.5end_60.fa > $BEDPE.species.no3Hetero.5end_60.information_content && \
#Rscript $PIPELINE_DIRECTORY/bin/draw_information_content.R  $BEDPE.species.no3Hetero.5end_60.information_content $BEDPE.species.no3Hetero.5end_60.information_content

rm -rf $BEDPE.reads.5end_60 $BEDPE.reads.5end_60.fa  $BEDPE.reads.5end_60.information_content $BEDPE.species.5end_60 $BEDPE.species.5end_60.fa $BEDPE.species.5end_60.information_content $BEDPE.species.no3Hetero.5end_60 $BEDPE.species.no3Hetero.5end_60.fa $BEDPE.species.no3Hetero.5end_60.information_content

