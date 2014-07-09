#!/bin/bash

export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis

#to process the output from /home/wangw1/git/smallRNA_analysis/g1Ut10A/pp8_q2_ww1_zscore_sep_0624_full.pl
#input file: eg. FLAGSiwiIP_KNOWNTE_AS.FLAGSiwiIP_KNOWNTE_AS.10.prefix.UA_VA.ppseq.txt
#trans	TGCTAGGGTTCGTGTTAGCAACGTCGT	nscaf2564,274886,+	CGTGTTAGCAACG	0.000514668039114771	AACCCTAGCAAGAGTCGTGCTTCGCAGA	scaffold14279,667,-	CGTGTTAGCATCG	0.00607287449392713	1111111111011

INDIR=/home/wangw1/uava/baseFraction/
OUTDIR=/home/wangw1/uava/baseFraction/unpaired
[ ! -f ${OUTDIR} ] && mkdir ${OUTDIR}

declare -a PIWI=("Ago3IPuniq" "AubIPuniq")
declare -a GROUPGT=("w1" "nosAgo3CDrescue" "aubvasAgo3CDrescue")

for gt in "${GROUPGT[@]}"
do 
	
	#total Ago3IPuniq
	cp /home/wangw1/data/projects/cd/ipsmRNA/jia_pipeline_results/Phil.SRA.Ago3IPuniq.${gt}.unox.ovary.inserts/Phil.SRA.Ago3IPuniq.${gt}.unox.ovary.inserts.xkxh.match2_all.out.23-29nt.uniq.reads.gz ${OUTDIR}
	gunzip ${OUTDIR}/Phil.SRA.Ago3IPuniq.${gt}.unox.ovary.inserts.xkxh.match2_all.out.23-29nt.uniq.reads.gz
	
	#total AubIPuniq
	cp /home/wangw1/data/projects/cd/ipsmRNA/jia_pipeline_results/Phil.SRA.AubIPuniq.${gt}.unox.ovary.inserts/Phil.SRA.AubIPuniq.${gt}.unox.ovary.inserts.xkxh.match2_all.out.23-29nt.uniq.reads.gz ${OUTDIR}
	gunzip ${OUTDIR}/Phil.SRA.AubIPuniq.${gt}.unox.ovary.inserts.xkxh.match2_all.out.23-29nt.uniq.reads.gz
	
	#cat cis and trans guide: Aub:Ago3 PP
	cat ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.cis.guide ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.trans.guide |cut -f1 |sort -u >${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.guide

	#cat cis and trans target: Aub:Ago3 PP
	cat ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.cis.target ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.trans.target |cut -f1 |sort -u >${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.target
	
	#unpaired guide: Aub
	exmatchBycol.pl ${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.guide ${OUTDIR}/Phil.SRA.AubIPuniq.${gt}.unox.ovary.inserts.xkxh.match2_all.out.23-29nt.uniq.reads 1 >${OUTDIR}/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.unpaired.guide

	#unpaired target: Ago3
	exmatchBycol.pl ${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.target ${OUTDIR}/Phil.SRA.Ago3IPuniq.${gt}.unox.ovary.inserts.xkxh.match2_all.out.23-29nt.uniq.reads 1 >${OUTDIR}/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.unpaired.target
	
				
	#cat cis and trans guide:Ago3:Aub PP
	cat ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.cis.guide ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.trans.guide |cut -f1 |sort -u >${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.guide

	#cat cis and trans target:Ago3:Aub PP
	cat ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.cis.target ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.trans.target |cut -f1 |sort -u >${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.target
	
	#unpaired guide: Ago3
	exmatchBycol.pl ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.guide ${OUTDIR}/Phil.SRA.Ago3IPuniq.${gt}.unox.ovary.inserts.xkxh.match2_all.out.23-29nt.uniq.reads 1 >${OUTDIR}/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.unpaired.guide	

	#unpaired guide: Aub
	exmatchBycol.pl ${INDIR}/${gt}_uniq_unox_AubIPSRA_Ago3IPSRA_prefix16_v26/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.target ${OUTDIR}/Phil.SRA.AubIPuniq.${gt}.unox.ovary.inserts.xkxh.match2_all.out.23-29nt.uniq.reads 1 >${OUTDIR}/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.unpaired.target	
																																																																																															
	echo "${PIPELINE_DIRECTORY}/Utils/base_fraction.pl -i ${OUTDIR}/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.unpaired.guide -o ${OUTDIR} -p 1 -r 2 -l 23" >>${OUTDIR}/pararun
	echo "${PIPELINE_DIRECTORY}/Utils/base_fraction.pl -i ${OUTDIR}/AubIPuniq_${gt}_unox_.Ago3IPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.unpaired.target -o ${OUTDIR} -p 1 -r 2 -l 23" >>${OUTDIR}/pararun
	echo "${PIPELINE_DIRECTORY}/Utils/base_fraction.pl -i ${OUTDIR}/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.unpaired.guide -o ${OUTDIR} -p 1 -r 2 -l 23" >>${OUTDIR}/pararun
	echo "${PIPELINE_DIRECTORY}/Utils/base_fraction.pl -i ${OUTDIR}/Ago3IPuniq_${gt}_unox_.AubIPuniq_${gt}_unox.16.prefix.UA_VA.ppseq.unpaired.target -o ${OUTDIR} -p 1 -r 2 -l 23" >>${OUTDIR}/pararun
done