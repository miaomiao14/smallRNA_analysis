#! /usr/bin/env bash
RC_EXT=200
GENOME_FA=/project/umw_phil_zamore/common/pipelines/piPipes/common/dm3/dm3.fa
CHROM=/project/umw_phil_zamore/common/pipelines/piPipes/common/dm3/dm3.ChromInfo.txt
CPU=5
PIPEDIR=/home/bh80w/Pipelines/piPipes/bin
outname=${2##*/} ##piRNAs
inname=${1##*/} ##degradome index
outdir=$3
echo -e ">deg" > ${1}r1.RC.ext${RC_EXT}.unique.fa && \
#${PIPEDIR}/
bedtools_piPipes slop -b $RC_EXT -i $1 -g $CHROM | \
awk -v total_len=$((2*RC_EXT+1)) 'BEGIN{OFS="\t"}{ if ($3-$2==total_len) { $6=($6=="+"?"-":"+"); print $0} }' | \
#${PIPEDIR}/
bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -tab -s -fo /dev/stdout | \
cut -f2 | \
gsort -u --temporary-directory=$PWD --parallel=$CPU >> ${1}r1.RC.ext${RC_EXT}.unique.fa && \
bowtie2-build ${1}r1.RC.ext${RC_EXT}.unique.fa ${1}r1.RC.ext${RC_EXT}.unique && \
awk -v piRNA_bot=23 -v piRNA_top=29 '{if (!printed[$7] && $3-$2 >= piRNA_bot && $3-$2 <= piRNA_top) {print ">"$7"_"$4"\n"$7; printed[$7]=1}}' $2 | \
bowtie2 \
	-N 0 -L 16 --gbar 16 --no-1mm-upfront -D 1 \
	-x ${1}r1.RC.ext${RC_EXT}.unique \
	-f -p $CPU \
	-U -  \
	1> ${outdir}/${outname%.bed*}.map2deg${inname}.sam \
	2> ${outdir}/${outname%.bed*}.map2deg${inname}.log && \
samtools view -bS ${outdir}/${outname%.bed*}.map2deg${inname}.sam | \
bedtools bamtobed -i stdin | \
awk -v len=$((2*RC_EXT+1)) '{if($6=="+") { a[$2%len]++; b[($3-1)%len]++} else {c[$2%len]++;d[($3-1)%len]++}}END{for (i=0;i<len;++i) printf "%d\t%d\t%d\t%d\t%d\n", i+1, a[i], b[i], -c[i], -d[i]}' \
	> ${outdir}/${outname%.bed*}.map2deg${inname}.species


awk -v piRNA_bot=23 -v piRNA_top=29 '{if (!printed[$7] && $3-$2 >= piRNA_bot && $3-$2 <= piRNA_top) { for (k=1; k<=$4; ++k) print ">"$7"_"$4"\n"$7; printed[$7]=1}}' $2 | \
bowtie2 \
	-N 0 -L 16 --gbar 16 --no-1mm-upfront -D 1 \
	-x ${1}r1.RC.ext${RC_EXT}.unique \
	-f -p $CPU \
	-U -  \
	1> ${outdir}/${outname%.bed*}.map2deg${inname}.sam \
	2> ${outdir}/${outname%.bed*}.map2deg${inname}.log && \
samtools view -bS ${outdir}/${outname%.bed*}.map2deg${inname}.sam | \
bedtools bamtobed -i stdin | \
awk -v len=$((2*RC_EXT+1)) '{if($6=="+") { a[$2%len]++; b[($3-1)%len]++} else {c[$2%len]++;d[($3-1)%len]++}}END{for (i=0;i<len;++i) printf "%d\t%d\t%d\t%d\t%d\n", i+1, a[i], b[i], -c[i], -d[i]}' \
	> ${outdir}/${outname%.bed*}.map2deg${inname}.reads
