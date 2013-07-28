#!/bin/bash
#input is the uniq.reads
#seq count
INDIR=$1
fileEND=$2
SPECIESLOG=${INDIR}/species.stat
READSLOG=${INDIR}/readsspecies.stat
paraFile=${INDIR}/species.reads.stat.${RANDOM}.paraFile
for i in `ls *.${fileEND}`
do 
echo -e "wc -l ${i} |awk 'BEGIN{OFS="\t"}{print $2,$1}' >> $SPECIESLOG " >> ${paraFile}
echo -e "awk -v filename=${i} 'BEGIN{OFS="\t"}{a+=\$2}END{print filename,a}' ${i} >> $READSLOG" >> ${paraFile}
done
if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands
fi
