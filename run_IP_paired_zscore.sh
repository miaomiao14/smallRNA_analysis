declare -a GROUPGT=("ago3cdmut_ox" "ago3wtmut_ox" "ago3cdwt_ox" "ago3cdmut_unox" "ago3wtmut_unox" "ago3cdwt_unox" \
"aubcdmut_ox" "aubwtmut_ox" "aubcdwt_ox" "aubcdmut_unox" "aubwtmut_unox" "aubcdwt_unox" \
"aubmuthet_ox" "aubhetmut_ox" "MTD_shaub_shGFP" "qinago3muthet_ox" "qinago3muthet_unox" \
"AubIP_aubcdwt_ox" "AubIP_aubcdwt_unox" \
)
declare -a ago3cdmut_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a ago3wtmut_ox=("Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a ago3cdwt_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts")

declare -a ago3cdmut_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a ago3wtmut_unox=("Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a ago3cdwt_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts")
#!/bin/bash


################
# Major Config #
################
# pipeline version
export smallRNA_downstream_analysis=1.0.0
# this pipeline is still in debug mode
export DEBUG=1 
# pipeline address: if you copy all the files to another directory, this is the place to change; under this directory sits two directories, bin and common. bin stores all the binary executables and common stores all the information of each ORGANISM.
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
# set PATH to be aware of pipeline/bin; this bin directory needs to be searched first
export PATH=${PIPELINE_DIRECTORY}/:$PATH

INDIR=$1 #this is the folder store all pipeline results outmost folders
OUT=${INDIR}/transposon_piRNA
LOG=${OUT}/log

declare -a GROUPGT=("AubIP_ago3cdwt_ox" "AubIP_ago3cdwt_unox" \
"AubIP_ago3cdw1_ox" "AubIP_ago3cdw1_unox" \
)

declare -a AubIP_ago3cdwt_ox=("Phil.AubIP.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.AubIP.aubvasAgo3WTrescue.ox.ovary.inserts")
declare -a AubIP_ago3cdwt_unox=("Phil.AubIP.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.AubIP.aubvasAgo3WTrescue.unox.ovary.inserts")

declare -a AubIP_ago3cdw1_ox=("Phil.AubIP.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.AubIP.w1.ox.ovary.inserts")
declare -a AubIP_ago3cdw1_unox=("Phil.AubIP.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.AubIP.w1.unox.ovary.inserts")

echo -e "`date` "+$ISO_8601"\tDraw paired transposon piRNA zscore scatterplot..." >> $LOG
OUTDIR10=${INDIR}/transposon_piRNA/paired_zscore_scatterplot
[ ! -d $OUTDIR10 ] && mkdir -p ${OUTDIR10}
[ ! -f ${OUT}/.status.${STEP}.transposon_piRNA.pairedZscore ] && \
paraFile=${OUTDIR10}/${RANDOM}.pairedZscore.para && \

for g in "${GROUPGT[@]}"
do
	SUBGROUP="$g[@]"
	[ -f ${OUTDIR10}/${g}.FB.zscore ] && rm ${OUTDIR10}/${g}.FB.zscore
	[ -f ${OUTDIR10}/${g}.total.zscore ] && rm ${OUTDIR10}/${g}.total.zscore
	count=1
	for t in ${!SUBGROUP}
	do
		cat ${INDIR}/pp6_FB/${t}/${t}.FB.pp6.out| awk -v gt=$t -v rank=$count '{OFS="\t"}{print gt,$1,$2,$3,$4,rank}' >> ${OUTDIR10}/${g}.FB.zscore
		count=$(($count+1))
		totalZscore=`cat ${INDIR}/pp6_FB/${t}/${t}.total.pp6.out|cut -f2`
		TZ=`printf "%0.2f" $totalZscore`
		echo -e ${t}"\t"${TZ} >>${OUTDIR10}/${g}.total.zscore
	done
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_zscore_FB_scatterplot ${OUTDIR10}/${g}.FB.zscore ${OUTDIR10}/${g}.total.zscore $g $OUTDIR10 " >>${paraFile}
done
[ $? == 0 ] && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands && \
	touch ${OUT}/.status.${STEP}.transposon_piRNA.pairedZscore
STEP=$((STEP+1))
