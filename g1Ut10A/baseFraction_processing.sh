#!/bin/sh

export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis

#to process the output from /home/wangw1/git/smallRNA_analysis/g1Ut10A/pp8_q2_ww1_zscore_sep_0624_full.pl
#input file: eg. FLAGSiwiIP_KNOWNTE_AS.FLAGSiwiIP_KNOWNTE_AS.10.prefix.UA_VA.ppseq.txt
#trans	TGCTAGGGTTCGTGTTAGCAACGTCGT	nscaf2564,274886,+	CGTGTTAGCAACG	0.000514668039114771	AACCCTAGCAAGAGTCGTGCTTCGCAGA	scaffold14279,667,-	CGTGTTAGCATCG	0.00607287449392713	1111111111011

INDIR=$1
OUTDIR=$2
[ ! -f ${OUTDIR} ] && mkdir -p ${OUTDIR}

declare -a GROUPGT=("cis" "trans")

for i in ${INDIR}/*.UA_VA.ppseq.txt
do
	filename=${i##*/}
	for g in "${GROUPGT[@]}"
	do
		cat $i |grep ${g} |cut -f2,3,4 |sort -u |cut -f1,3 >${OUTDIR}/${filename%.txt}.${g}.guide
		${PIPELINE_DIRECTORY}/Utils/base_fraction.pl -i ${OUTDIR}/${filename%.txt}.${g}.guide -o ${OUTDIR} -p 1 -r 2 -l 23
		
		cat $i |grep ${g} |cut -f5,6,7 |sort -u |cut -f1,3 >${OUTDIR}/${filename%.txt}.${g}.target
		${PIPELINE_DIRECTORY}/Utils/base_fraction.pl -i ${OUTDIR}/${filename%.txt}.${g}.target -o ${OUTDIR} -p 1 -r 2 -l 23
	done
done