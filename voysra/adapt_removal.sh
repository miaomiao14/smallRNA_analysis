export PATH=${wd}/${project_folder}/:$HOME/bin:$PATH
project_folder=$1
wd=/diag/home/netashawang/scratch/
mkdir -p ${wd}/${project_folder}/adapter_trim/tooshort
for INPUT in ${wd}/${project_folder}/raw/*.fastq.gz
do
IFILE=${INPUT##*/}
ODIR=${INPUT%/*}
OFILE=${ODIR%/*}
echo "cutadapt -a TGGAATTCTCGGGTGCCAAGG -e 0.1 --no-indels -O 8 --trimmed-only -m 15 -M 50  $INPUT --too-short-output=${wd}/${project_folder}/adapter_trim/tooshort/${IFILE/fastq.gz/}tooshort.fastq  --trim-n -q 20,20 >${wd}/${project_folder}/adapter_trim/${IFILE/fastq.gz/}trim.fastq && gzip ${wd}/${project_folder}/adapter_trim/${IFILE/fastq.gz/}trim.fastq" >>${wd}/${project_folder}/adapt_removal.run.sh
done
chmod 755 ${wd}/${project_folder}/adapt_removal.run.sh
~/src/HpcGridRunner-1.0.2/hpc_cmds_GridRunner.pl -G ~/sge.conf -c ${wd}/${project_folder}/adapt_removal.run.sh
