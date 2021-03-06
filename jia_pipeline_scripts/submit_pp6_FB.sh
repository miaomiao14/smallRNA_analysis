#!/bin/bash

export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis

###This is for heterotypic analysis with pp6
## input: file1 file2 dir

INDIR=$1 #this is the folder store all pipeline results outmost folders

#declare -a TARGETS=("Julius.SRA.MTD_shaub.unox.ovary.inserts" "Julius.SRA.MTD_shGFP.unox.ovary.inserts" \
#"Phil.AubIP.AubCDrescue.ox.ovary.inserts" "Phil.AubIP.AubCDrescue.unox.ovary.inserts" \
#"Phil.AubIP.AubWTrescue.ox.ovary.inserts" "Phil.AubIP.AubWTrescue.unox.ovary.inserts" \
#"Phil.SRA.aubHets.ox.ovary.inserts" "Phil.SRA.aubMutant.ox.ovary.inserts" \
#"Phil.SRA.QinAgo3Hets.ox.ovary.inserts" "Phil.SRA.QinAgo3Hets.unox.ovary.inserts" \
#"Phil.SRA.QinAgo3Muts.ox.ovary.inserts" "Phil.SRA.QinAgo3Muts.unox.ovary.inserts")

#t="Phil.Ago3IPuniq.aubMutant.unox.ovary.inserts"
declare -a TARGETS=("Phil.SRA.aubMutsrep2.unox.ovary.inserts" "Phil.SRA.ago3Mutsrep2.unox.ovary.inserts" "Phil.SRA.AubWTrescuerep2.unox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescuerep2.ox.ovary.inserts" "Phil.SRA.AubCDrescuerep2.unox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescuerep2.ox.ovary.inserts" "Phil.SRA.ago3Mutsrep2.ox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescuerep2.unox.ovary.inserts" "Phil.SRA.w1CJ.ox.ovary.inserts" "Phil.SRA.AubWTrescuerep2.ox.ovary.inserts" "Phil.SRA.aubMutsrep2.ox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescuerep2.unox.ovary.inserts" "Phil.SRA.AubCDrescuerep2.ox.ovary.inserts" "Phil.SRA.ago3MutsAubMuts.ox.ovary.inserts" "Phil.SRA.ago3MutsAubMuts.unox.ovary.inserts")
[ ! -d ${INDIR}/pp6_FB ] && mkdir -p ${INDIR}/pp6_FB
for t in ${TARGETS[@]}
do
for i in `ls ${INDIR}/${t}/*.xkxh.transposon.mapper2.gz`
#for i in `ls ${INDIR}/*.inserts/*.xkxh.transposon.mapper2.gz`
do 
	#ln -s $i ${OUTDIR}
	FILE=${i##*/}
	FILENAME=${FILE%.gz}
	insertsname=`basename $FILE .xkxh.transposon.mapper2.gz`
	#inserts=${FILE%%.inserts.*}
	#inserts=${inserts}.inserts
	#nfnnc=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	#nfdep=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	OUTDIR=${INDIR}/pp6_FB/${insertsname}
	mkdir -p ${OUTDIR}
	SGE=${INDIR}/pp6_FB/${insertsname}.pp6_FB.sge

	script=${PIPELINE_DIRECTORY}/pp6_T_ww_len.pl
	
#echo "spliting transposon family..."
	echo "#!/bin/sh

#$ -pe single 8
#$ -V
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash
##$ -M weiwanghhq@gmail.com
#$ -m e 

	export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis

	[ ! -f $OUTDIR/${insertsname}.total.pp6.out ] && $script $i $i 1 ${OUTDIR} >$OUTDIR/${insertsname}.total.pp6.out &&
	${PIPELINE_DIRECTORY}/FB.pl $i ${OUTDIR} &&
	paraFile=${OUTDIR}/${insertsname}.pp6.para
	for j in \`ls -1 ${OUTDIR}/${FILENAME}.*\`
	do
	T1=\${j##*mapper2.}
	
	"echo -ne  \" T1=\${j##*mapper2.} \&\& \" \>\> \${paraFile}"
	#"echo -ne \" \\\`$script \${j} \${j} 1 ${OUTDIR} \\\$T1 \\\` \>\> $OUTDIR/${insertsname}.FB.\${T1}.pp6.out \&\& \" \>\> \${paraFile}"
	#"echo -e \" \\\`$script \${j} \${j} 1 ${OUTDIR} \\\$T1 \\\` \>\> $OUTDIR/${insertsname}.FB.pp6.temp  \" \>\> \${paraFile}"
	"echo -ne \" $script \${j} \${j} 1 ${OUTDIR} \\\$T1 \>\> $OUTDIR/${insertsname}.FB.\${T1}.pp6.out \&\& \" \>\> \${paraFile}"
	"echo -e \" $script \${j} \${j} 1 ${OUTDIR} \\\$T1 \>\> $OUTDIR/${insertsname}.FB.pp6.temp  \" \>\> \${paraFile}"

	done
	if [[ ! -f \${paraFile}.completed ]] || [[ -f \$paraFile.failed_commands ]]
	then
	
		ParaFly -c \$paraFile -CPU 24 -failed_cmds \$paraFile.failed_commands
	fi
	awk '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${OUTDIR}/${insertsname}.FB.pp6.temp > ${OUTDIR}/${insertsname}.FB.pp6.out
"> $SGE
	
	qsub $SGE
	#sleep 5s

done
done