#! /usr/bin/env bash
#input: inserts1 inserts2 <inserts3>

export PIPELINEDIR=/home/lees2/pipeline:/home/xuj1/pipeline

INPUTDIR=/home/wangw1/isilon_temp/ipsmRNA/trimmed
OUTPUTDIR=/home/wangw1/isilon_temp/ipsmRNA/trimmed_uniqBound
[ ! -d ${OUTPUTDIR} ] && mkdir -p ${OUTPUTDIR}
declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue")
declare -a OX=("ox" "unox")
echo -e "genotype\tox\tsharedSpecies\tsharedReadsAgo3IP\tsharedReadsAubIP\tuniqSpeciesAgo3IP\t" >> ${OUTPUTDIR}/stat.LOG
echo -e "uniqSpeciesAubIP\tuniqReadsAgo3IP\tuniqReadsAubIP\n" >> ${OUTPUTDIR}/${g}.${o}.log
for g in "${GT[@]}"
do
	for o in "${OX[@]}"
	do
	cd ${OUTPUTDIR}	
	A=Phil.SRA.Ago3IP.${g}.${o}.ovary.trimmed
	B=Phil.SRA.AubIP.${g}.${o}.ovary.trimmed
	ln -s ${INPUTDIR}/${A}
	ln -s ${INPUTDIR}/${B}
 
	/home/wangw1/bin/inserts_file_methods.py -a $A -b $B -i
	eval sharedSpecies=`wc -l ${A}.sharedA|cut -d" " -f1`
	eval sharedReadsAgo3IP=`sumcol ${A}.sharedA`
	eval sharedReadsAubIP=`sumcol ${B}.sharedB`
	
	/home/wangw1/bin/inserts_file_methods.py -a $A -b $B -u 
	
	uniqSpeciesAgo3IP=`wc -l ${A}.uniqA`
	uniqSpeciesAubIP=`wc -l ${B}.uniqB`
	uniqReadsAgo3IP=`sumcol ${A}.uniqA` 
	uniqReadsAubIP=`sumcol ${B}.uniqB`

	echo -e "${g}\t${o}\t${sharedSpecies}\t${sharedReadsAgo3IP}\t${sharedReadsAubIP}\t${uniqSpeciesAgo3IP}\t${uniqSpeciesAubIP}\t${uniqReadsAgo3IP}\t${uniqReadsAubIP}\n" >> ${OUTPUTDIR}/stat.LOG

	done
done