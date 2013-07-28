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
echo -ne "wc -l ${i} >> $SPECIESLOG " >> ${paraFile}
echo -e "awk -v filename=${i} '{a+=\$2}END{print a,filename}' ${i} >> $READSLOG" >> ${paraFile}
done
if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands
fi
