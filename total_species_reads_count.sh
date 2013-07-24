#!/bin/bash

INDIR=$1
fileEND=$2
LOG=${INDIR}/species.reads.stat
paraFile=${INDIR}/species.reads.stat.${RANDOM}.paraFile
for i in `ls *.${fileEND}`
do 
echo -ne "wc -l ${i} >> $LOG &&" >> ${paraFile}
echo -e "awk -v filename=${i} '{a+=\$2}END{print a,filename}' ${i} >> $LOG" >> ${paraFile}
done
if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands
fi
