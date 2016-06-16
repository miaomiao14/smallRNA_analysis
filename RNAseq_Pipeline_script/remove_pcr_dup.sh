#!/bin/sh

FILE1=$1
FILE2=$2
FILE=${1##*/}
#prefix=`basename $FILE .[left|right].fq`
prefix=`echo -e "${FILE1}\n${FILE2}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && prefix=${prefix%.*}
prefix=`echo $prefix | sed -e s/[^a-z0-9]*$//g`
#####Prepare fastq files by removing qual scores and placing the sequence as a 2nd column to the IDs.
[ ! -f ${prefix}.read1.noqual.txt ] && grep /1$ -A 1 $FILE1 | sed '/--/d' | awk '/[0-9]$/{printf $0" ";next;}1' | sed 's#/1##g' > ${prefix}.read1.noqual.txt
[ ! -f ${prefix}.read2.noqual.txt ] && grep /2$ -A 1 $FILE2 | sed '/--/d' | awk '/[0-9]$/{printf $0" ";next;}1' | sed 's#/2##g' > ${prefix}.read2.noqual.txt

#####Merge sequences and remove duplicates.
awk '{p=$2; getline<f} !A[$2=p $2]++' f=${prefix}.read2.noqual.txt ${prefix}.read1.noqual.txt | tr " " "\t" > ${prefix}.reads_rmdup.txt && \
#####Prepare read IDs to be used for searching.
cut -f1 ${prefix}.reads_rmdup.txt > ${prefix}.reads_rmdup_IDs.txt && \
sed 's#$#/1#g' ${prefix}.reads_rmdup_IDs.txt > ${prefix}.read1_rmdup_IDs.txt && \
sed 's#$#/2#g' ${prefix}.reads_rmdup_IDs.txt > ${prefix}.read2_rmdup_IDs.txt && \

#####Get non-duplicate fastq reads.
awk 'NR==FNR{A[$1]=$1;next} $1 in A{print;getline;print;getline;print;getline;print}' FS="\t" ${prefix}.read1_rmdup_IDs.txt FS="\t" OFS="\t" $FILE1 > ${prefix}.read1_rmdup.fastq && \
awk 'NR==FNR{A[$1]=$1;next} $1 in A{print;getline;print;getline;print;getline;print}' FS="\t" ${prefix}.read2_rmdup_IDs.txt FS="\t" OFS="\t" $FILE2 > ${prefix}.read2_rmdup.fastq && \
rm ${prefix}.read1.noqual.txt && \
rm ${prefix}.read2.noqual.txt && \
rm ${prefix}.reads_rmdup.txt && \
rm ${prefix}.reads_rmdup_IDs.txt && \
rm ${prefix}.read1_rmdup_IDs.txt 
rm ${prefix}.read2_rmdup_IDs.txt
