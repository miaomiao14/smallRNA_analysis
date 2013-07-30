#!/bin/bash

export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis

###This is for heterotypic analysis with pp6
## input: file1 file2 dir

INDIR=$1 #this is the folder store all pipeline results outmost folders
[ ! -d ${INDIR}/pp6_FB ] && mkdir -p ${INDIR}/pp6_FB

for i in `ls ${INDIR}/*.inserts/*.xkxh.transposon.mapper2.gz`
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

	script=${PIPELINE_DIRECTORY}/pp6_T_ww.pl
	
#echo "spliting transposon family..."
	echo "#!/bin/sh

#$ -pe single 24
#$ -V
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash
##$ -M weiwanghhq@gmail.com
#$ -m e 

export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis

	${PIPELINE_DIRECTORY}/FB.pl \$i \${OUTDIR} &&
	paraFile=${OUTDIR}/${insertsname}.pp6.para
	for j in \`ls -1 \${OUTDIR}/\${FILENAME}.*\`
	do
	T1=\${j##*mapper2.}
	
	"echo -ne  \" T1=\\\${j##*mapper2.} \&\& \" \>\> ${paraFile}"
	"echo -e \" \\\`$script \\\${j} \\\${j} 1 ${OUTDIR} \\\$T1 \\\` \>\> $OUTDIR/${insertsname}.FB.\${T1}.pp6.out  \" \>\> \${paraFile}"
	"echo -e \" \\\`$script \\\${j} \\\${j} 1 ${OUTDIR} \\\$T1 \\\` \>\> $OUTDIR/${insertsname}.FB.pp6.temp  \" \>\> \${paraFile}"

	done
	if [[ ! -f \${paraFile}.completed ]] || [[ -f \$paraFile.failed_commands ]]
	then
	
		ParaFly -c \$paraFile -CPU 24 -failed_cmds \$paraFile.failed_commands
	fi
"> $SGE
	
	#qsub $SGE
	#sleep 5s
#awk '{OFS="\t"}{print $0}' ${OUTDIR}/${insertsname}.FB.pp6.temp > ${OUTDIR}/${insertsname}.FB.pp6
done
 