#!/bin/sh
#A=$1
#B=$2
F1=${1##*/}
F2=${2##*/}
OUTDIR=$3
n=$4

for i in {1..25}
do
echo "#!/bin/sh
#$ -V
#$ -pe single 8
#$ -o \$HOME/sge_jobs_output/sge_job.\$JOB_ID.out -j y
#$ -S /bin/bash
#$ -m e 

export PIPELINEDIR=/home/lees2/pipeline:/home/xuj1/pipeline:/home/xuj1/pipeline_mmu
export PATH=\$PATH:\$PIPELINEDIR
export PATH=/share/bin/R/bin:/share/bin/R/share:\$PATH
export LD_LIBRARY_PATH=/share/bin/R/lib64:\$LD_LIBRARY_PATH

cd $OUTDIR
mkdir ${i}_basepair
cd $OUTDIR/${i}_basepair
ln -s $1 $F1
N1=${F1%*.xkxh.norm.bed}
if [ ! -e $F2 ]
then
ln -s $2 $F2
N2=${F2%*.xkxh.norm.bed}
else
N2=\$N1
fi
/home/wangw1/bin/pp8_q3_ww3_revised_06172012_percentage.pl $F1 $F2 $n $i >> $OUTDIR/\${N1}_\${N2}_extensive_complementarity.log 
">submit_pp8_${i}bp_overlap.sh
qsub submit_pp8_${i}bp_overlap.sh
done
