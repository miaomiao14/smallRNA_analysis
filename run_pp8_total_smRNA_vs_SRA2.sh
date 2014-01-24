#!/bin/bash -x
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
script=${PIPELINE_DIRECTORY}/pp8_ww_smRNA_vs_DEG.pl

OUT0=/home/wangw1/isilon_temp/smRNA/pp8_smRNA_vs_smRNA_total_01242014


#c=$6 #cpu
declare -a GROUPGT=("/scratch/hanb/cd/smallRNA_pipeline_output/Brennecke.SRA.WT.unox.OSC/intersect_piRNA_length/Brennecke.SRA.WT.unox.OSC.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Brennecke.SRA.zucMut.unox.ovary/intersect_piRNA_length/Brennecke.SRA.zucMut.unox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.Ago3IPsds.aubMuts.unox.ovary/intersect_piRNA_length/Phil.SRA.Ago3IPsds.aubMuts.unox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.Ago3IPsds.w1.unox.ovary/intersect_piRNA_length/Phil.SRA.Ago3IPsds.w1.unox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.ago3MutsAubMuts.ox.ovary/intersect_piRNA_length/Phil.SRA.ago3MutsAubMuts.ox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.ago3Mutsrep2.ox.ovary/intersect_piRNA_length/Phil.SRA.ago3Mutsrep2.ox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.AubCDrescuerep2.ox.ovary/intersect_piRNA_length/Phil.SRA.AubCDrescuerep2.ox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.aubMutsrep2.ox.ovary/intersect_piRNA_length/Phil.SRA.aubMutsrep2.ox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.aubvasAgo3CDrescuerep2.ox.ovary/intersect_piRNA_length/Phil.SRA.aubvasAgo3CDrescuerep2.ox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.aubvasAgo3WTrescuerep2.ox.ovary/intersect_piRNA_length/Phil.SRA.aubvasAgo3WTrescuerep2.ox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.AubWTrescuerep2.ox.ovary/intersect_piRNA_length/Phil.SRA.AubWTrescuerep2.ox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2" \
"/scratch/hanb/cd/smallRNA_pipeline_output/Phil.SRA.w1.ox.ovary/intersect_piRNA_length/Phil.SRA.w1.ox.ovary.trimmed.x_rRNA.x_hairpin.x_dmel_virus.dm3v0.uniq.L23mer.x_rpmk_rtRNA.bed2")

declare -a FEATURE=("all")

STEP=1

[ ! -d ${OUT0} ] && mkdir -p ${OUT0}
#step 1
#OUT=${OUT0}

if [ ! -f ${OUT0}/.status.${STEP}.SRA_DEG.pp8 ] 
then
for gt in "${GROUPGT[@]}"
do
	for f in "${FEATURE[@]}"
	do
		filename=${gt##*/}
		gf=${filename%.trimmed*}
		g=${gf#*SRA.}
		OUTDIR=${OUT0}/${g}_${f}
		[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
		
		LOG=${OUTDIR}/${g}.log	
		smmapper2=${OUT0}/${gf}.inserts.xkxh.all.norm.bed.23-29 #share the SRA norm.bed files
		[ ! -f $smmapper2 ] && \
		awk 'BEGIN{OFS="\t"}{if($3-$2<30) print $1,$2,$3,$6,$7,$4,$5}' $gt >${smmapper2} && gzip ${smmapper2}.gz	
	#total Ping-Pong
		[ ! -s ${OUT0}/${g}.total.pp8.out ] && submitsge 8 ${g} ${OUT0} "$script ${smmapper2}.gz ${smmapper2}.gz 1 ${OUTDIR} >${OUT0}/${g}.total.pp8.out" 
	done

done
fi
[ $? == 0 ] && \
touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8
STEP=$((STEP+1))

declare -a GROUPGT=c("zucMut" "w1" "AubWTrescuerep2" "aubvasAgo3WTrescuerep2" "aubvasAgo3CDrescuerep2" "ago3Mutsrep2" "aubMutsrep2" "AubCDrescuerep2" "ago3MutsAubMuts")
#generate master table for ppscore
masterOUT=${OUT0}/masterpp8score
[ ! -d ${masterOUT} ] && mkdir -p ${masterOUT}
if [ ! -f ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master ] 
then
	for f in "${FEATURE[@]}"
	do
		[ -f ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.txt ] && rm ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.txt

		
		[ -f ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.txt ] && rm ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.txt

		for g in "${GROUPGT[@]}"
		do
			#filename=${gt##*/}
			#gf=${filename%.trimmed*}
			#g=${gf#*SRA.}
			OUTDIR=${OUT0}/${g}_${f}
			cut -f1,2 ${OUTDIR}/${g}_SRA_all.${g}_SRA_all.pp|awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}' >> ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.txt	
			cut -f1,3 ${OUTDIR}/${g}_SRA_all.${g}_SRA_all.pp|awk -v gt=$g 'BEGIN{OFS="\t"}{print gt,$1,$2}' >> ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.txt
	
		done
		
	
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.txt ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.mastertable.txt
		
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.txt ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.mastertable.txt
	done
fi

[ $? == 0 ] && \
	touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master
STEP=$((STEP+1))