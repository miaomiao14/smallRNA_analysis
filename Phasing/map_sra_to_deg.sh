#! /usr/bin/env bash
RC_EXT=200
GENOME_FA=/home/wangw1/data/common/dm3.TAS.fa
CHROM=/home/wangw1/data/common/dm3.chromInfo
CPU=5
PIPEDIR=~/git/smallRNA_analysis/Phasing
echo -e ">deg" > ${1}r1.RC.ext${RC_EXT}.unique.fa && \
${PIPEDIR}/bedtools_piPipes slop -b $RC_EXT -i $1 -g $CHROM | \
awk -v total_len=$((2*RC_EXT+1)) 'BEGIN{OFS="\t"}{ if ($3-$2==total_len) { $6=($6=="+"?"-":"+"); print $0} }' | \
${PIPEDIR}/bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -tab -s -fo /dev/stdout | \
cut -f2 | \
gsort -u --temporary-directory=$PWD --parallel=$CPU >> ${1}r1.RC.ext${RC_EXT}.unique.fa && \
bowtie2-build ${1}r1.RC.ext${RC_EXT}.unique.fa ${1}r1.RC.ext${RC_EXT}.unique && \
awk -v piRNA_bot=23 -v piRNA_top=29 '{if (!printed[$7] && $3-$2 >= piRNA_bot && $3-$2 <= piRNA_top) {print ">"$7"_"$4"\n"$7; printed[$7]=1}}' $2 | \
bowtie2 \
	-x ${1}r1.RC.ext${RC_EXT}.unique \
	-f -p $CPU \
	-U -  \
	1> ${2%.bed*}.map2deg${1}.sam \
	2> ${2%.bed*}.map2deg${1}.log && \
samtools view -bS ${2%.bed*}.map2deg${1}.sam | \
/home/wangw1/wangwj2bin/bedtools bamtobed -i stdin | \
awk -v len=$((2*RC_EXT+1)) '{if($6=="+") { a[$2%len]++; b[($3-1)%len]++} else {c[$2%len]++;d[($3-1)%len]++}}END{for (i=0;i<len;++i) printf "%d\t%d\t%d\t%d\t%d\n", i+1, a[i], b[i], -c[i], -d[i]}' \
	> ${2%.bed*}.map2deg${1}.specie