#! /usr/bin/env bash

#input file is xkxh.norm.bed or xkxh.transposon.mapper2
#give input directory, 
#binsize is for the bigWigSummary
INDIR=$1 #this is the folder store all pipeline results outmost folders
BINSIZE=$2 
[ ! -d ${INDIR}/circos ] && mkdir -p ${INDIR}/circos

for i in `ls ${INDIR}/*.inserts/*.xkxh.transposon.mapper2.gz`
do 
	#ln -s $i ${OUTDIR}
	FILE=${i##*/}
	insertsname=`basename $FILE .xkxh.transposon.mapper2.gz`
	inserts=${FILE%%.inserts.*}
	inserts=${inserts}.inserts
	nfnnc=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	nfdep=`cat ${INDIR}/${inserts}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	OUTDIR=${INDIR}/circos
	mkdir -p ${OUTDIR}
	SGE=${INDIR}/circos/${insertsname}.circos.bin${BINSIZE}.sge
	
###generate sge file
	echo "#!/bin/sh
#$ -V
#$ -cwd
#$ -pe single 8
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -l mem_free=31G
#$ -S /bin/bash
#$ -m e
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
export PATH=${PIPELINE_DIRECTORY}/:$PATH
mkdir \$HOME/scratch/jobid_\$JOB_ID
export DIR=\$HOME/scratch/jobid_\$JOB_ID
start=\$SGE_TASK_FIRST
end=\$SGE_TASK_LAST

FILETYPE=mapper2
CHRSIZE=/home/wangw1/pipeline/common/dm3.chrom.sizes

declare -a NORMFACTOR=(${nfnnc} ${nfdep})
declare -a NORMFACTORTYPE=("nnc" "seqDep")
count=0
for NF in \${NORMFACTOR[@]}
do 
\${PIPELINE_DIRECTORY}/normbedmapper2circos.pl ${i} \$FILETYPE \$CHRSIZE \$NF \${NORMFACTORTYPE[\$count]} \$DIR
paraFile=\${DIR}/\${RANDOM}.bed2bw.para

	for j in \`ls \${DIR}/*circos.bed\`
	do 
		"echo -ne \" bedSort \$j \$j.sort \&\& \"  \>\>\${paraFile}"
		"echo -ne \" bedItemOverlapCountWithScore dm3 \$j.sort chromSize=\$CHRSIZE stdin \>\${j}.bedGraph \&\& \"  \>\>\${paraFile}"
		"echo -e \" bedGraphToBigWig \${j}.bedGraph \$CHRSIZE \${j}.bw \"  \>\>\${paraFile}"
	done
	if [[ ! -f \${paraFile}.completed ]] || [[ -f \$paraFile.failed_commands ]]
		then
	
		ParaFly -c \$paraFile -CPU 8 -failed_cmds \$paraFile.failed_commands
	fi
	
	for j in \`ls \${DIR}/*circos.bed\`
	do
		a=(\$(cat \$CHRSIZE))
		b=$BINSIZE
		paraFile=\${DIR}/\${RANDOM}.bw2bwsummary.para
		for k in \$(seq 0 2 \$(( \${#a[@]} - 1)))
		do 
			tlen=\${a[\$((\$k+1))]}
			n=\$((\$tlen/\$b))
			nbin=\$((\$n+1))
			
			"echo -e \" bigWigSummary \${j}.bw \${a[\$k]} 0 \$tlen \$nbin -type=mean \> \${j}.\${a[\$k]}.mean.\${b}.txt \"  \>\>\${paraFile}"
			"echo -e \" bigWigSummary \${j}.bw \${a[\$k]} 0 \$tlen \$nbin -type=max \> \${j}.\${a[\$k]}.max.\${b}.txt \" \>\> \${paraFile}"
		done
		if [[ ! -f \${paraFile}.completed ]] || [[ -f \$paraFile.failed_commands ]]
		then
	
		ParaFly -c \$paraFile -CPU 8 -failed_cmds \$paraFile.failed_commands
		fi
	
	done

paraFile=\${DIR}/\${RANDOM}.bwSummary2bincount.para
"echo -e \" \${PIPELINE_DIRECTORY}/bedscorefile2circostrack.pl \$DIR \${NORMFACTORTYPE[\$count]} sense mean \$b ${i} \$FILETYPE \$DIR \" \>\> \${paraFile} "
"echo -e \" \${PIPELINE_DIRECTORY}/bedscorefile2circostrack.pl \$DIR \${NORMFACTORTYPE[\$count]} antisense mean \$b ${i} \$FILETYPE \$DIR \" \>\> \${paraFile} "
"echo -e \" \${PIPELINE_DIRECTORY}/bedscorefile2circostrack.pl \$DIR \${NORMFACTORTYPE[\$count]} sense max \$b ${i} \$FILETYPE \$DIR \" \>\> \${paraFile} "
"echo -e \" \${PIPELINE_DIRECTORY}/bedscorefile2circostrack.pl \$DIR \${NORMFACTORTYPE[\$count]} antisense max \$b ${i} \$FILETYPE \$DIR \" \>\> \${paraFile} "

if [[ ! -f \${paraFile}.completed ]] || [[ -f \$paraFile.failed_commands ]]
then
	
	ParaFly -c \$paraFile -CPU 8 -failed_cmds \$paraFile.failed_commands
fi
count=\$((\$count+1))
done

rm -f \${DIR}\*.bin.txt 
rm -f \${DIR}\*.\${BINSIZE}.txt
rm -f \${DIR}\*.bed
rm -f \${DIR}\*.bedGraph
rm -f \${DIR}\*.sort
rm -f \${DIR}\*.bw
[ -d ${OUTDIR}/${insertsname}.$BINSIZE ] && rm -rf ${OUTDIR}/${insertsname}.$BINSIZE
mv \$HOME/scratch/jobid_\$JOB_ID ${OUTDIR}/${insertsname}.$BINSIZE

"> $SGE
qsub $SGE
done


