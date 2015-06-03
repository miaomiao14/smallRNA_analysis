#! /usr/bin/env bash
#input: inserts1 inserts2 <inserts3>

export PIPELINEDIR=/home/lees2/pipeline:/home/xuj1/pipeline

INPUTDIR=/home/wangw1/isilon_temp/ipsmRNA/trimmed
OUTPUTDIR=/home/wangw1/isilon_temp/ipsmRNA/trimmed_uniqBoundsds
[ ! -d $OUTPUTDIR ] && mkdir $OUTPUTDIR

STEP=1

[ ! -d ${OUTPUTDIR} ] && mkdir -p ${OUTPUTDIR}
declare -a GT=("w1" "AubCDrescue" "AubWTrescue" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue")
#declare -a OX=("ox" "unox")
#declare -a GT=("ago3Hets" "aubHets" "qinHets" "nosAgo3CDrescue" "nosAgo3WTrescue")
declare -a OX=("unox")
declare -a PRO=("Ago3IPsds" "AubIP")
echo -e "genotype\tox\tsharedSpecies\tsharedReadsAgo3IPsds\tsharedReadsAubIP\tuniqSpeciesAgo3IPsds\tuniqSpeciesAubIP\tuniqReadsAgo3IPsds\tuniqReadsAubIP\n" >> ${OUTPUTDIR}/stat.log
parafile=${OUTPUTDIR}/para.uniq.bound
[ -s ${parafile} ] && rm ${parafile}
[ ! -f ${OUTPUTDIR}/.status.${STEP}.IPuniqBound ] && \
for g in "${GT[@]}"
do
	for o in "${OX[@]}"
	do
	cd ${OUTPUTDIR}	
	echo -ne " A=Phil.SRA.${PRO[0]}.${g}.${o}.ovary.trimmed && " >>${parafile}
	echo -ne " B=Phil.SRA.${PRO[1]}.${g}.${o}.ovary.trimmed && " >>${parafile}
	echo -ne " ANEW=Phil.SRA.${PRO[0]}uniq.${g}.${o}.ovary.trimmed && " >>${parafile}
	echo -ne " BNEW=Phil.SRA.${PRO[1]}uniq.${g}.${o}.ovary.trimmed && " >>${parafile}
	echo -ne " ln -s ${INPUTDIR}/\${A} ${OUTPUTDIR} && " >>${parafile}
	echo -ne " ln -s ${INPUTDIR}/\${B} ${OUTPUTDIR} && " >>${parafile}
 
	echo -ne " /home/wangw1/bin/inserts_file_methods.py -a \$A -b \$B -i && " >>${parafile}
	echo -ne " sharedSpecies=\`wc -l \${A}.sharedA|cut -d\" \" -f1\` && " >>${parafile}
	echo -ne " sharedReadsAgo3IPsds=\`sumcol \${A}.sharedA 2\` && " >>${parafile}
	echo -ne " sharedReadsAubIP=\`sumcol \${B}.sharedB 2\` && " >>${parafile}
	
	#echo -ne " A=Phil.SRA.Ago3IPsds.${g}.${o}.ovary.trimmed && " >>${parafile}
	#echo -ne " B=Phil.SRA.AubIP.${g}.${o}.ovary.trimmed && " >>${parafile}
	echo -ne "/home/wangw1/bin/inserts_file_methods.py -a \$A -b \$B -u && " >>${parafile}
	
	echo -ne " uniqSpeciesAgo3IPsds=\`wc -l \${A}.uniqA|cut -d\" \" -f1\` && " >>${parafile}
	echo -ne " uniqSpeciesAubIP=\`wc -l \${B}.uniqB|cut -d\" \" -f1\` && " >>${parafile}
	echo -ne " uniqReadsAgo3IPsds=\`sumcol \${A}.uniqA 2\` && " >>${parafile} 
	echo -ne " uniqReadsAubIP=\`sumcol \${B}.uniqB 2\` && rm \${A} && rm \${B} && mv \${A}.uniqA \${ANEW} && mv \${B}.uniqB \${BNEW} &&  " >>${parafile}

	echo -e " echo -e \"${g}\\\t${o}\\\t\${sharedSpecies}\\\t\${sharedReadsAgo3IPsds}\\\t\${sharedReadsAubIP}\\\t\${uniqSpeciesAgo3IPsds}\\\t\${uniqSpeciesAubIP}\\\t\${uniqReadsAgo3IPsds}\\\t\${uniqReadsAubIP}\" >> ${OUTPUTDIR}/stat.log " >>${parafile}

	done
done
[ $? == 0 ] && \
if [[ ! -f ${parafile}.completed ]] || [[ -f $parafile.failed_commands ]]
then
	CPUN=`wc -l $parafile |cut -f1 -d" "` && \
	ParaFly -c $parafile -CPU $CPUN -failed_cmds $parafile.failed_commands
fi
[ $? == 0 ] && \
touch ${OUTPUTDIR}/.status.${STEP}.IPuniqBound
STEP=$((STEP+1))