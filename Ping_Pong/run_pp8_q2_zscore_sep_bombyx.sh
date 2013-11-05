#!/bin/sh
#A=$1
#B=$2
F1=${1##*/}
F2=${2##*/}
OUTDIR=$3
n=$4
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
ln -s $1 $F1
N1=${F1%*.xkxh.norm.bed.gz}
if [ ! -e $F2 ]
then
ln -s $2 $F2
N2=${F2%*.xkxh.norm.bed.gz}
else
N2=\$N1
fi
/home/wangw1/bin/pp8_q2_ww1_zscore_sep.pl $F1 $F2 $n  > $OUTDIR/\${N1}_\${N2}_pp8_q2_UA_VA.log 
">submit_UA_VA.sh
qsub submit_UA_VA.sh

totalBED=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.bed2
submitsge 24 TOTALBEDPP6 ~/isilon_temp/BmN4/Yuki.SRA.TOTAL.DMSO.ox.BmN4cell.inserts/pp6_total "awk '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$4/\$5,\$6}' ${totalBED} >${totalBED%.bed2}.normalized.bed && /home/wangw1/git/smRNA_analysis/pp6_ww_bed.pl -i ${totalBED%.bed2}.normalized.bed -o ~/isilon_temp/BmN4/Yuki.SRA.TOTAL.DMSO.ox.BmN4cell.inserts/pp6_total -f bedscore "

AGO3BED=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.bed2
SIWIBED=/home/wangw1/isilon_temp/BmN4/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.inserts/Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell.bmv2v0.all.all.xrRNA.xtRNA.xh.bed2
submitsge 24 AGO3SIWIPP6 ~/isilon_temp/BmN4/pp6_TOTAL "awk '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$4/\$5,\$6}' ${AGO3BED} >${AGO3BED%.bed2}.normalized.bed && awk '{OFS=\"\\t\"}{print \$1,\$2,\$3,\$4,\$4/\$5,\$6}' ${SIWIBED} >${SIWIBED%.bed2}.normalized.bed && \
/home/wangw1/git/smRNA_analysis/pp6_ww_bed.pl -i ${AGO3BED%.bed2}.normalized.bed -j ${SIWIBED%.bed2}.normalized.bed -o ~/isilon_temp/BmN4/pp6_TOTAL -f bedscore"


declare -a GT=("Yuki.SRA.FLAGBmAgo3IP.DMSO.ox.BmN4cell" "Yuki.SRA.FLAGSiwiIP.DMSO.ox.BmN4cell" "Yuki.SRA.TOTAL.DMSO.ox.BmN4cell")
for t in ${GT[@]}
do
	BED=/home/wangw1/isilon_temp/BmN4/${t}.inserts/${t}.bmv2v0.all.all.xrRNA.xtRNA.xh.bed2
	awk '{OFS="\t"}{print $1,$2,$3,$4,$4/$5,$6}' ${BED} >${BED%.bed2}.normalized.bed
	/home/wangw1/bin/submitsge 24 ${t} $OUTDIR "/home/wangw1/git/smRNA_analysis/pp6_ww_bed.pl -i ${totalBED%.bed2}.normalized.bed -o $OUTDIR -f bedscore"
done

#converting to xkxh.norm.bed format, run pp8_q2_ww1_zscore_sep.pl

#zip the mapper2 format, run pp6
#Ago3IP S: SiwiIP AS, Ago3IP AS: SiwiIP S, Ago3IP total: SiwiIP total	