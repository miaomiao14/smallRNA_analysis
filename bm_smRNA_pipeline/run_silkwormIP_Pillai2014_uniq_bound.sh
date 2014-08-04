#! /usr/bin/env bash
#input: inserts1 inserts2 <inserts3>

export PIPELINEDIR=/home/lees2/pipeline:/home/xuj1/pipeline

INPUTDIR=/home/wangw1/data/bmn4/BmN4/Xiol2014
OUTPUTDIR=/home/wangw1/data/bmn4/BmN4/Xiol2014/uniqBound
[ ! -d $OUTPUTDIR ] && mkdir $OUTPUTDIR

STEP=1

[ ! -d ${OUTPUTDIR} ] && mkdir -p ${OUTPUTDIR}
declare -a GT=("wt")
#declare -a OX=("ox" "unox")
#declare -a GT=("ago3Hets" "aubHets" "qinHets" "nosAgo3CDrescue" "nosAgo3WTrescue")
declare -a OX=("unox")
declare -a PRO=("Ago3IP" "SiwiIP")
echo -e "genotype\tox\tsharedSpecies\tsharedReadsAgo3IP\tsharedReadsSiwiIP\tuniqSpeciesAgo3IP\tuniqSpeciesSiwiIP\tuniqReadsAgo3IP\tuniqReadsSiwiIP\n" >> ${OUTPUTDIR}/stat.log
parafile=${OUTPUTDIR}/para.uniq.bound
[ -s ${parafile} ] && rm ${parafile}
[ ! -f ${OUTPUTDIR}/.status.${STEP}.IPuniqBound ] && \
for g in "${GT[@]}"
do
	for o in "${OX[@]}"
	do
	cd ${OUTPUTDIR}	
	echo -ne " A=Pillai2014.SRA.${PRO[0]}.${g}.${o}.BmN4cell.trimmed && " >>${parafile}
	echo -ne " B=Pillai2014.SRA.${PRO[1]}.${g}.${o}.BmN4cell.trimmed && " >>${parafile}
	echo -ne " ANEW=Pillai2014.SRA.${PRO[0]}uniq.${g}.${o}.BmN4cell.trimmed && " >>${parafile}
	echo -ne " BNEW=Pillai2014.SRA.${PRO[1]}uniq.${g}.${o}.BmN4cell.trimmed && " >>${parafile}
	echo -ne " ln -s ${INPUTDIR}/\${A} ${OUTPUTDIR} && " >>${parafile}
	echo -ne " ln -s ${INPUTDIR}/\${B} ${OUTPUTDIR} && " >>${parafile}
 
	echo -ne " /home/wangw1/bin/inserts_file_methods.py -a \$A -b \$B -i && " >>${parafile}
	echo -ne " sharedSpecies=\`wc -l \${A}.sharedA|cut -d\" \" -f1\` && " >>${parafile}
	echo -ne " sharedReadsAgo3IP=\`sumcol \${A}.sharedA 2\` && " >>${parafile}
	echo -ne " sharedReadsSiwiIP=\`sumcol \${B}.sharedB 2\` && " >>${parafile}
	
	#echo -ne " A=Pillai2014.SRA.Ago3IP.${g}.${o}.BmN4cell.trimmed && " >>${parafile}
	#echo -ne " B=Pillai2014.SRA.SiwiIP.${g}.${o}.BmN4cell.trimmed && " >>${parafile}
	echo -ne "/home/wangw1/bin/inserts_file_methods.py -a \$A -b \$B -u && " >>${parafile}
	
	echo -ne " uniqSpeciesAgo3IP=\`wc -l \${A}.uniqA|cut -d\" \" -f1\` && " >>${parafile}
	echo -ne " uniqSpeciesSiwiIP=\`wc -l \${B}.uniqB|cut -d\" \" -f1\` && " >>${parafile}
	echo -ne " uniqReadsAgo3IP=\`sumcol \${A}.uniqA 2\` && " >>${parafile} 
	echo -ne " uniqReadsSiwiIP=\`sumcol \${B}.uniqB 2\` && rm \${A} && rm \${B} && mv \${A}.uniqA \${ANEW} && mv \${B}.uniqB \${BNEW} &&  " >>${parafile}

	echo -e " echo -e \"${g}\\\t${o}\\\t\${sharedSpecies}\\\t\${sharedReadsAgo3IP}\\\t\${sharedReadsSiwiIP}\\\t\${uniqSpeciesAgo3IP}\\\t\${uniqSpeciesSiwiIP}\\\t\${uniqReadsAgo3IP}\\\t\${uniqReadsSiwiIP}\" >> ${OUTPUTDIR}/stat.log " >>${parafile}

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