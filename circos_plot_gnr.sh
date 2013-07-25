#! /usr/bin/env bash

#input file is xkxh.norm.bed or xkxh.transposon.mapper2
#give input directory, 
#binsize is for the bigWigSummary
INDIR=$1 #this is the folder store all pipeline results outmost folders
BINSIZE=$2 
[ ! -d ${INDIR}/circos ] && mkdir -p ${INDIR}/circos

for i in `ls ${INDIR}/*.inserts/*uniqmap.xkxh.transposon.mapper2.gz`
do 
	ln -s $i ${OUTDIR}
	FILE=\${i##*/}
	insertsname=`basename $FILE .xkxh.transposon.mapper2.gz`
	nfnnc=`cat ${INDIR}/${insertsname}/output/${insertsname}_stats_table_reads|tail -1|cut -f4`
	nfdep=`cat ${INDIR}/${insertsname}/output/${insertsname}_stats_table_reads|tail -1|cut -f2`
	OUTDIR=${INDIR}/circos/${insertsname}
	SGE=${INDIR}/circos/${insertsname}.circos.sge
	
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

declar -a NORMFACTOR=(${nfnnc} ${nfdep})

for NF in ${NORMFACTOR[@]}
do \
${PIPELINE_DIRECTORY}/normbedmapper2circos.pl \${i} \$FILETYPE \$CHRSIZE $NF \$DIR
	for j in \`ls \${DIR}/*circos.bed\`
	do \
		bedSort \$j \$j.sort
		bedItemOverlapCountWithScore dm3 \$j.sort chromSize=\$CHRSIZE stdin >\${j}.bedGraph
		bedGraphToBigWig \${j}.bedGraph \$CHRSIZE \${j}.bw
		a=(\$(cat \$CHRSIZE))
		b=$BINSIZE
		for k in \$(seq 0 2 \$(( \${#a[@]} - 1)))
		do \
			tlen=\${a[\$((\$k+1))]}
			n=\$((\$tlen/\$b))
			nbin=\$((\$n+1))
			bigWigSummary \${j}.bw \${a[\$k]} 0 \${a[\$((\$k+1))]} \$nbin -type=mean >\${j}.\${a[\$k]}.mean.bin.txt
			bigWigSummary \${j}.bw \${a[\$k]} 0 \${a[\$((\$k+1))]} \$nbin -type=max >\${j}.\${a[\$k]}.max.bin.txt
		done
	
	done
done

#match2all_scorefile2circostrack.pl \$DIR mean \$bin \$file
#match2all_scorefile2circostrack.pl \$DIR max \$bin \$file

#\`rm -f *.bin.txt\` 
#\`rm -f *.bed\`
#\`rm -f *.bedGraph\`
#\`rm -f *.sort\`
#\`rm -f *.bw\`
#for k in \`ls \${DIR}/*.\$bin.txt\`
#do
#echo \$k
#cut -f4 \$k |sort -r -n |head -1
#done
mv \$HOME/scratch/jobid_\$JOB_ID ${OUTDIR}/${insertsname}
"> $SGE
qsub $SGE
done


