#! /usr/bin/env bash
#input file is match2all.out
#binsize is for the bigWigSummary
Dir=$1
File=$2
binsize=$3
#mapper2circos.pl ${DIR}/$file $chrsize
#for i in `ls ${DIR}/*sense.bed`
#do
#F=${i##*/}
SGE=${File}.circos.sge
echo "#!/bin/sh
#$ -V
#$ -cwd
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -l mem_free=10G
#$ -S /bin/bash
#$ -m e
# pipeline address: if you copy all the files to another directory, this is the place to change; under this directory sits two directories, bin and common. bin stores all the binary executables and common stores all the information of each ORGANISM.
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
# set PATH to be aware of pipeline/bin; this bin directory needs to be searched first
export PATH=${PIPELINE_DIRECTORY}/:$PATH
export RUNDIR=/home/wangw1/bin
mkdir \$HOME/scratch/jobid_\$JOB_ID
start=\$SGE_TASK_FIRST
end=\$SGE_TASK_LAST

DIR=$Dir
file=$File
chrsize=/home/wangw1/pipeline/common/dm3.chrom.sizes
NF=$normfile
lenrangeselector \${DIR}/\$file 23 29 > \${DIR}/\${file}.23-29
/home/wangw1/bin/normbed2circos.pl \${file}.23-29 \$chrsize $NF
echo \"file format done!\"
for i in \`ls \${DIR}/*circos.bed\`
do
F=\${i##*/}
bedSort \$i \$F.sort
bedItemOverlapCountWithScore dm3 \$F.sort chromSize=\$chrsize stdin >\${F}.bedGraph
bedGraphToBigWig \${F}.bedGraph \$chrsize \${F}.bw
a=(\$(cat \$chrsize))
#echo " Number of elements in array is \$\(\( $\{#a[\@]\} \)\)"
bin=$binsize
for j in \$(seq 0 2 \$(( \${#a[@]} - 1)))
do
tlen=\${a[\$((\$j+1))]}
n=\$((\$tlen/\$bin))
nbin=\$((\$n+1))
bigWigSummary \${F}.bw \${a[\$j]} 0 \${a[\$((\$j+1))]} \$nbin -type=mean >\${F}.\${a[\$j]}.mean.bin.txt
bigWigSummary \${F}.bw \${a[\$j]} 0 \${a[\$((\$j+1))]} \$nbin -type=max >\${F}.\${a[\$j]}.max.bin.txt
done
done

match2all_scorefile2circostrack.pl \$DIR mean \$bin \$file
match2all_scorefile2circostrack.pl \$DIR max \$bin \$file

\`rm -f *.bin.txt\` 
\`rm -f *.bed\`
\`rm -f *.bedGraph\`
\`rm -f *.sort\`
\`rm -f *.23-29\`
\`rm -f *.bw\`
for k in \`ls \${DIR}/*.\$bin.txt\`
do
echo \$k
cut -f4 \$k |sort -r -n |head -1
done
"> $SGE
qsub $SGE