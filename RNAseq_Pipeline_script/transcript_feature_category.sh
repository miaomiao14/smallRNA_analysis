#!/bin/bash 

#test
#This script is to assign piRNAs to different transcript features
#all exons
#First exon(5«UTR) (]
#Last exon(3«UTR)  [)
#middle exon []
#all introns ][
#Exon-Intron junctions
#Intron-Exon junctions
#Exon-Exon junctions

export PIPELINE_DIRECTORY=/diag/home/netashawang/git/RNAseq_Pipeline

#The input is the gtf transcript file from cufflinks
GTF=$1
OUTDIR=$2
CHRSIZE=$3 #name length
FASTA=$4 #genome fasta file
FILE=${1##*/}
prefix=${FILE%.gtf}
LOG=${OUTDIR}/${prefix}.log
gtfToGenePred $GTF ${OUTDIR}/${prefix}.gp 
genePredToBed ${OUTDIR}/${prefix}.gp ${OUTDIR}/${prefix}.bed
rm ${OUTDIR}/${prefix}.gp
bedtools sort -i ${OUTDIR}/${prefix}.bed > ${OUTDIR}/${prefix}.sorted.bed && \
rm ${OUTDIR}/${prefix}.bed 
bedClip ${OUTDIR}/${prefix}.sorted.bed $CHRSIZE ${OUTDIR}/${prefix}.sorted.clipped.bed && \
rm ${OUTDIR}/${prefix}.sorted.bed
bedToGenePred ${OUTDIR}/${prefix}.sorted.clipped.bed  ${OUTDIR}/${prefix}.sorted.clipped.gp
bedToTxEdges ${OUTDIR}/${prefix}.sorted.clipped.bed ${OUTDIR}/${prefix}.sorted.clipped.txEdges

paraFile=${OUTDIR}/${RANDOM}.para
echo -e "awk 'BEGIN{OFS=\"\\\t\"}{if( \$7==\"(\" && \$9==\"]\" ) print \$1,\$2,\$3,\$4,\$5,\$6 }'  ${OUTDIR}/${prefix}.sorted.clipped.txEdges |sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse |awk 'BEGIN{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$5,"0",\$4}' > ${OUTDIR}/${prefix}.sorted.clipped.firstexon.bed" >$paraFile
echo -e "awk 'BEGIN{OFS=\"\\\t\"}{if( \$7==\"[\" && \$9==\")\" ) print \$1,\$2,\$3,\$4,\$5,\$6 }'  ${OUTDIR}/${prefix}.sorted.clipped.txEdges |sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse |awk 'BEGIN{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$5,"0",\$4}'> ${OUTDIR}/${prefix}.sorted.clipped.lastexon.bed" >>$paraFile
echo -e "awk 'BEGIN{OFS=\"\\\t\"}{if( \$7==\"[\" && \$9==\"]\" ) print \$1,\$2,\$3,\$4,\$5,\$6 }'  ${OUTDIR}/${prefix}.sorted.clipped.txEdges |sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse |awk 'BEGIN{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$5,"0",\$4}'> ${OUTDIR}/${prefix}.sorted.clipped.middleexon.bed" >>$paraFile
echo -e "awk 'BEGIN{OFS=\"\\\t\"}{if( \$7==\"]\" && \$9==\"[\" ) print \$1,\$2,\$3,\$4,\$5,\$6 }'  ${OUTDIR}/${prefix}.sorted.clipped.txEdges |sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse|awk 'BEGIN{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$5,"0",\$4}' > ${OUTDIR}/${prefix}.sorted.clipped.intron.bed" >>$paraFile
echo -e "awk 'BEGIN{OFS=\"\\\t\"}{if( \$8==\"exon\" ) print \$1,\$2,\$3,\$4,\$5,\$6 }'  ${OUTDIR}/${prefix}.sorted.clipped.txEdges |sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse|awk 'BEGIN{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$5,"0",\$4}' > ${OUTDIR}/${prefix}.sorted.clipped.allexon.bed" >>$paraFile
echo -e "awk 'BEGIN{OFS=\"\\\t\"}{if( \$7==\"(\" && \$9==\")\" ) print \$1,\$2,\$3,\$4,\$5,\$6 }'  ${OUTDIR}/${prefix}.sorted.clipped.txEdges |sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse |awk 'BEGIN{OFS=\"\\\t\"}{print \$1,\$2,\$3,\$5,"0",\$4}'> ${OUTDIR}/${prefix}.sorted.clipped.singleexon.bed" >>$paraFile
CPU=`wc -l $paraFile|cut -f1`
ParaFly -c $paraFile -CPU $CPU -failed_cmds $paraFile.failed_commands
FLANKING=35
#Exon-Intron Junction
#use exon ends position
[ ! -f ${OUTDIR}/${prefix}.sorted.clipped.ExonIntronSplice.bed ] && awk -v flank=$FLANKING 'BEGIN{OFS="\t"}{if($8>1){num=split ($10,b,",");for(i=1;i<num;i++) {start=b[i]-flank;end=b[i]+flank;if(start>0) {print $2,start,end,$1,"0",$3"\n"}}}}' ${OUTDIR}/${prefix}.sorted.clipped.gp|sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,"0",$4}' >${OUTDIR}/${prefix}.sorted.clipped.ExonIntronSplice.bed && \
bedtools slop -g $CHRSIZE -b 0 -i ${OUTDIR}/${prefix}.sorted.clipped.ExonIntronSplice.bed >${OUTDIR}/${prefix}.sorted.clipped.ExonIntronSplice.slopped.bed
#Intron-Exon Junction
#use exon starts position
[ ! -f ${OUTDIR}/${prefix}.sorted.clipped.IntronExonSplice.bed ] && awk -v flank=$FLANKING 'BEGIN{OFS="\t"}{if($8>1){num=split ($9,a,",");for(i=2;i<=num;i++) {start=a[i]-flank;end=a[i]+flank;if(start>0) {print $2,start,end,$1,"0",$3"\n"}}}}' ${OUTDIR}/${prefix}.sorted.clipped.gp|sort -k1,1 -k2,2n -k3,3n -k6,6 | groupBy -g 1,2,3,6 -c 4, -o collapse |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,"0",$4}' >${OUTDIR}/${prefix}.sorted.clipped.IntronExonSplice.bed && \
bedtools slop -g $CHRSIZE -b 0 -i ${OUTDIR}/${prefix}.sorted.clipped.IntronExonSplice.bed >${OUTDIR}/${prefix}.sorted.clipped.IntronExonSplice.slopped.bed
#Exon-Exon Junction
#use exon end position and exon start position
[ ! -f ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.bed ] && awk -v flank=$FLANKING 'BEGIN{OFS="\t"}{if($8>1){num=split ($9,a,",");ends=split ($10,b,",");for(i=1;i<num;i++) {start1=b[i]-flank;end1=b[i];start2=a[i];end2=a[i]+flank;if(start1>0) {print $2,start1,end1,$1"_EEJ"i,"0",$3"\n"$2,start2,end2,$1"_EEJ"i,"0",$3"\n"}}}}' ${OUTDIR}/${prefix}.sorted.clipped.gp >${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.bed && \
bedtools slop -g $CHRSIZE -b 0 -i ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.bed >${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.bed && \
rm ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.bed
[ ! -f ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.fa ] && fastaFromBed -tab -s -name -fi $FASTA -bed ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.bed -fo stdout | awk '{OFS="\n"}{seq1=$2;getline;seq2=$2;print ">"$1,seq1""seq2}' >${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.fa
[ ! -f ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.1.ebwt ] && bowtie-build ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped.fa ${OUTDIR}/${prefix}.sorted.clipped.ExonExonJun.slopped

#calculate the genome coverage of each annotation feature, print it to the lOG file
${PIPELINE_DIRECTORY}/feature_coverage.sh $CHRSIZE ${OUTDIR} ${prefix} $LOG

