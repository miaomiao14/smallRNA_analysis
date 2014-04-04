#! /usr/bin/env bash
#input: inserts1 inserts2 <inserts3>

export PIPELINEDIR=/home/lees2/pipeline:/home/xuj1/pipeline

INPUTDIR=/home/wangw1/data/mouse/
OUTPUTDIR=/home/wangw1/data/mouse/uniqBound
[ ! -d $OUTPUTDIR ] && mkdir $OUTPUTDIR

STEP=1

[ ! -d ${OUTPUTDIR} ] && mkdir -p ${OUTPUTDIR}
declare PI="Carroll"
declare -a GT=("Hets")
declare -a ORG=("mice")
declare -a OX=("unox")
declare -a PRO=("MiliIP" "Miwi2IP")
declare -a REP=("rep1" "rep2")
echo -e "genotype\tox\treplicate\tsharedSpecies\tsharedReadsMiliIP\tsharedReadsMiwi2IP\tuniqSpeciesMiliIP\tuniqSpeciesMiwi2IP\tuniqReadsMiliIP\tuniqReadsMiwi2IP\n" >> ${OUTPUTDIR}/stat.log
parafile=${OUTPUTDIR}/para.uniq.bound
[ -s ${parafile} ] && rm ${parafile}
[ ! -f ${OUTPUTDIR}/.status.${STEP}.IPuniqBound ] && \
for g in "${GT[@]}"
do
	for o in "${OX[@]}"
	do
		for r in "${REP[@]}"
		do
		cd ${OUTPUTDIR}	
		echo -ne " A=${PI}.SRA.${PRO[0]}.${g}.${o}.${ORG}.${r}.mm0.trimmed && " >>${parafile}
		echo -ne " B=${PI}.SRA.${PRO[1]}.${g}.${o}.${ORG}.${r}.mm0.trimmed && " >>${parafile}
		echo -ne " ANEW=${PI}.SRA.${PRO[0]}uniq.${g}.${o}.${ORG}.${r}.mm0.trimmed && " >>${parafile}
		echo -ne " BNEW=${PI}.SRA.${PRO[1]}uniq.${g}.${o}.${ORG}.${r}.mm0.trimmed && " >>${parafile}
		echo -ne " ln -s ${INPUTDIR}/\${A} ${OUTPUTDIR} && " >>${parafile}
		echo -ne " ln -s ${INPUTDIR}/\${B} ${OUTPUTDIR} && " >>${parafile}
	 
		echo -ne " /home/wangw1/bin/inserts_file_methods.py -a \$A -b \$B -i && " >>${parafile}
		echo -ne " sharedSpecies=\`wc -l \${A}.sharedA|cut -d\" \" -f1\` && " >>${parafile}
		echo -ne " sharedReadsMiliIP=\`sumcol \${A}.sharedA 2\` && " >>${parafile}
		echo -ne " sharedReadsMiwi2IP=\`sumcol \${B}.sharedB 2\` && " >>${parafile}
		
	#echo -ne " A=${PI}.SRA.MiliIP.${g}.${o}.${ORG}.trimmed && " >>${parafile}
	#echo -ne " B=${PI}.SRA.Miwi2IP.${g}.${o}.${ORG}.trimmed && " >>${parafile}
		echo -ne "/home/wangw1/bin/inserts_file_methods.py -a \$A -b \$B -u && " >>${parafile}
		
		echo -ne " uniqSpeciesMiliIP=\`wc -l \${A}.uniqA|cut -d\" \" -f1\` && " >>${parafile}
		echo -ne " uniqSpeciesMiwi2IP=\`wc -l \${B}.uniqB|cut -d\" \" -f1\` && " >>${parafile}
		echo -ne " uniqReadsMiliIP=\`sumcol \${A}.uniqA 2\` && " >>${parafile} 
		echo -ne " uniqReadsMiwi2IP=\`sumcol \${B}.uniqB 2\` && rm \${A} && rm \${B} && mv \${A}.uniqA \${ANEW} && mv \${B}.uniqB \${BNEW} &&  " >>${parafile}
	
		echo -e " echo -e \"${g}\\\t${r}\\\t${o}\\\t\${sharedSpecies}\\\t\${sharedReadsMiliIP}\\\t\${sharedReadsMiwi2IP}\\\t\${uniqSpeciesMiliIP}\\\t\${uniqSpeciesMiwi2IP}\\\t\${uniqReadsMiliIP}\\\t\${uniqReadsMiwi2IP}\" >> ${OUTPUTDIR}/stat.log " >>${parafile}
		done
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