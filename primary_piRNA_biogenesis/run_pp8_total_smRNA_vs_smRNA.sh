#!/bin/bash -x
export PIPELINE_DIRECTORY=/home/ww74w/git/smallRNA_analysis/
script=${PIPELINE_DIRECTORY}/primary_piRNA_biogenesis/pp8_ww_smRNA_vs_DEG_ghpcc.pl

OUT0=/home/bh80w/scratch/SRA/pp8

declare -a FEATURE=("all")

STEP=1

[ ! -d ${OUT0} ] && mkdir -p ${OUT0}
#step 1
#OUT=${OUT0}
INPUT=/home/bh80w/scratch/SRA/all_x_rpmk_MASK_bed2
if [ ! -f ${OUT0}/.status.${STEP}.SRA_SRA.pp8 ] 
then
for i in ${INPUT}/*.bed2
	
do
	gt=${i%.trimmed*}
	for f in "${FEATURE[@]}"
	do
		filename=${gt##*/}
		gf=${filename%.trim*}
		#g=${gf#*SRA.}
		g=$gf
		OUTDIR=${OUT0}/${g}_${f}
		[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
		
		LOG=${OUTDIR}/${g}.log	
		smmapper2=${OUT0}/${gf}.inserts.xkxh.all.norm.bed.23-29 #share the SRA norm.bed files
		[ ! -f $smmapper2 ] && \
			
	#total Ping-Pong
[ ! -s ${OUT0}/${g}.total.pp8.out ] && /home/ww74w/smallRNApipeline/bin/hm.bsub.wrapper 1 10:00 short ${g} ${OUT0} "awk 'BEGIN{OFS=\"\\t\"}{if(\$3-\$2<30) {start=\$2+1;print \$1,start,\$3,\$6,\$7,\$4,\$5}}' $i >${smmapper2} && gzip ${smmapper2} && $script ${smmapper2}.gz ${smmapper2}.gz 1 ${OUTDIR} >${OUT0}/${g}.total.pp8.out" 
	done

done
fi
[ $? == 0 ] && \
touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8
STEP=$((STEP+1))

#declare -a GROUPGT=("zucMut" "w1" "AubWTrescuerep2" "aubvasAgo3WTrescuerep2" "aubvasAgo3CDrescuerep2" "ago3Mutsrep2" "aubMutsrep2" "AubCDrescuerep2" "ago3MutsAubMuts")
#generate master table for ppscore

masterOUT=${OUT0}/masterpp8score
[ ! -d ${masterOUT} ] && mkdir -p ${masterOUT}
touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master
if [ ! -f ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master ] 
then
	for f in "${FEATURE[@]}"
	do
		[ -f ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.txt ] && rm ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.txt

		
		[ -f ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.txt ] && rm ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.txt

		for gt in "${GROUPGT[@]}"
		do
			filename=${gt##*/}
			gf=${filename%.trimmed*}
			gn=${gf#*SRA.}
			g=${gn%%.[ox|unox]*}
			
			OUTDIR=${OUT0}/${gn}_${f}
			cut -f1,2 ${OUTDIR}/${g}_SRA_all.${g}_SRA_all.pp|awk -v gt=$gn 'BEGIN{OFS="\t"}{print gt,$1,$2}' >> ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.txt	
			cut -f1,3 ${OUTDIR}/${g}_SRA_all.${g}_SRA_all.pp|awk -v gt=$gn 'BEGIN{OFS="\t"}{print gt,$1,$2}' >> ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.txt
	
		done
		
	
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.txt ${masterOUT}/SRA_all.SRA_all.nonnormalized.pp8score.mastertable.txt
		
	${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source cast_master_table ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.txt ${masterOUT}/SRA_all.SRA_all.normalized.pp8score.mastertable.txt
	done
fi

[ $? == 0 ] && \
	touch ${OUT0}/.status.${STEP}.SRA_DEG.pp8.master
STEP=$((STEP+1))