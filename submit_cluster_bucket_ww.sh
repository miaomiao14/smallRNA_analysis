#!/bin/bash

# USAGE: make_sge_cluster_bucket.sh <n> <input xkxh.norm.bed file1 with full path>  <stats table1> <input xkxh.norm.bed file2 with full path> <stats table2> <n> <option [empty: default Brennecke 142 ovary clusters;  custom defined cluster file]>
n=$1
F1=${2##*/}
S1=${F1%.xkxh.norm.bed}
#S1=${F1%.uniqmap}
F2=${4##*/}
S2=${F2%.xkxh.norm.bed}
#S2=${F2%.uniqmap}

#S= ${S1}_${S2}

#if [ $# -ne 7 ] 
#then 
#C='/home/xuj1/nearline/cluster/data/JB.allcluster'
#awk '{print $1,"(+)"}' $C | sed 's/ //g'  > $2/$S.cluster.all
#else
#C=$7
#awk '{print $6,"(+)"}' $C | sed 's/ //g'  > $2/$S.cluster.all
#fi

echo "#!/bin/sh
#$ -V
#$ -pe single 4
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash

export PIPELINEDIR=/home/xuj1/nearline/cluster/data:/home/xuj1/pipeline:/home/xuj1/pipeline/bowtie:/home/lees2/pipeline:
export OUTPUTDIR=$6
export PATH=\$PATH:\$PIPELINEDIR
export PATH=/share/bin/R/bin:/share/bin/R/share:\$PATH
export LD_LIBRARY_PATH=/share/bin/R/lib64:\$LD_LIBRARY_PATH
mkdir \$HOME/scratch/jobid_\$JOB_ID
#cp $1 $2/$S.cluster.all \$HOME/scratch/jobid_\$JOB_ID
cd \$HOME/scratch/jobid_\$JOB_ID
ln -s $2 $F1
ln -s $4 $F2
cp $3 .
cp $5 .
#/home/xuj1/nearline/cluster/cluster_bucket/cluster.convertseq.pl $S.cluster.all > $S.cluster.fa
#bowtie-build $S.cluster.fa $S
/home/wangw1/bin/cluster_bucket3_ww.pl /home/wangw1/pipeline/common/JB.allcluster.fa $n $F1 $3 $F2 $5
cp \$HOME/scratch/jobid_\$JOB_ID/${F1}_${F2}.pdf \$OUTPUTDIR 
#rm $OUTPUTDIR/$S.cluster.all 
"> submit_${F1}_${F2}_cluster.sge
qsub submit_${F1}_${F2}_cluster.sge

