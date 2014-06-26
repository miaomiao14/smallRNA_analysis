#!/bin/sh

export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis

#to process the output from /home/wangw1/git/smallRNA_analysis/g1Ut10A/pp8_q2_ww1_zscore_sep_0624_full.pl
#input file: eg. FLAGSiwiIP_KNOWNTE_AS.FLAGSiwiIP_KNOWNTE_AS.10.prefix.UA_VA.ppseq.txt
#trans	TGCTAGGGTTCGTGTTAGCAACGTCGT	nscaf2564,274886,+	CGTGTTAGCAACG	0.000514668039114771	AACCCTAGCAAGAGTCGTGCTTCGCAGA	scaffold14279,667,-	CGTGTTAGCATCG	0.00607287449392713	1111111111011

INDIR=/home/wangw1/data/projects/cd/ipsmRNA/jia_pipeline_results
OUTDIR=/home/wangw1/uava/baseFraction_cistrans/total
[ ! -f ${OUTDIR} ] && mkdir ${OUTDIR}

for gt in `find ${INDIR} -wholename "*uniq*inserts.xkxh.match2_all.out.23-29nt.uniq.reads.gz"`
do

	echo "${PIPELINE_DIRECTORY}/Utils/base_fraction.pl -i ${gt} -o ${OUTDIR} -p 1 -r 2 -l 23" >>${OUTDIR}/pararun

done