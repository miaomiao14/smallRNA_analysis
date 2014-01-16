#! /usr/bin/env bash
#input: inserts1 inserts2 <inserts3>

export PIPELINEDIR=/home/lees2/pipeline:/home/xuj1/pipeline

INPUTDIR=$1
OUTPUTDIR=$2
[ ! -d $OUTPUTDIR ] && mkdir $OUTPUTDIR

STEP=1

declare -a GT=("GT1" "GT2")
declare -a GT1=("nosGal4Cyo" "nosGal4UASPArmiWT" )
declare -a GT2=("nosGal4Cyo" "nosGal4UASPArmiK729A")
declare -a LEN=("30" "78B")
declare -a OX=("unox")
declare -a PRO=("FLAGIPnon" "FLAGIP")
echo -e "genotype\tox\tsharedSpecies\tsharedReadsFLAGIPnon\tsharedReadsFLAGIP\tuniqSpeciesFLAGIPnon\tuniqSpeciesFLAGIP\tuniqReadsFLAGIPnon\tuniqReadsFLAGIP\n" >> ${OUTPUTDIR}/stat.log
parafile=${OUTPUTDIR}/para.uniq.bound
[ -s ${parafile} ] && rm ${parafile}
[ ! -f ${OUTPUTDIR}/.status.${STEP}.IPuniqBound ] && \
for g in "${GT[@]}"
do
	eval "SUBGROUP=(\"\${${g}[@]}\")"
	for o in "${OX[@]}"
	do
		for l in "${LEN[@]}"
		do
			cd ${OUTPUTDIR}	
			echo -ne " A=Phil.SRA.${PRO[0]}.${SUBGROUP[0]}${l}.${o}.ovary.trimmed && " >>${parafile}
			echo -ne " B=Phil.SRA.${PRO[1]}.${SUBGROUP[1]}${l}.${o}.ovary.trimmed && " >>${parafile}
			echo -ne " ANEW=Phil.SRA.${PRO[0]}uniq.${SUBGROUP[0]}${l}.${SUBGROUP[1]}${l}.${o}.ovary.trimmed && " >>${parafile}
			echo -ne " BNEW=Phil.SRA.${PRO[1]}uniq.${SUBGROUP[1]}${l}.${SUBGROUP[0]}${l}.${o}.ovary.trimmed && " >>${parafile}
			echo -ne " ln -s ${INPUTDIR}/\${A} ${OUTPUTDIR} && " >>${parafile}
			echo -ne " ln -s ${INPUTDIR}/\${B} ${OUTPUTDIR} && " >>${parafile}
		 
			echo -ne " /home/wangw1/bin/inserts_file_methods.py -a \$A -b \$B -i && " >>${parafile}
			echo -ne " sharedSpecies=\`wc -l \${A}.sharedA|cut -d\" \" -f1\` && " >>${parafile}
			echo -ne " sharedReadsFLAGIPnon=\`sumcol \${A}.sharedA 2\` && " >>${parafile}
			echo -ne " sharedReadsFLAGIP=\`sumcol \${B}.sharedB 2\` && " >>${parafile}
			

			echo -ne "/home/wangw1/bin/inserts_file_methods.py -a \$A -b \$B -u && " >>${parafile}
			
			echo -ne " uniqSpeciesFLAGIPnon=\`wc -l \${A}.uniqA|cut -d\" \" -f1\` && " >>${parafile}
			echo -ne " uniqSpeciesFLAGIP=\`wc -l \${B}.uniqB|cut -d\" \" -f1\` && " >>${parafile}
			echo -ne " uniqReadsFLAGIPnon=\`sumcol \${A}.uniqA 2\` && " >>${parafile} 
			echo -ne " uniqReadsFLAGIP=\`sumcol \${B}.uniqB 2\` && rm \${A} && rm \${B} && mv \${A}.uniqA \${ANEW} && mv \${B}.uniqB \${BNEW} &&  " >>${parafile}
		
			echo -e " echo -e \"${g}\\\t${o}\\\t\${sharedSpecies}\\\t\${sharedReadsFLAGIPnon}\\\t\${sharedReadsFLAGIP}\\\t\${uniqSpeciesFLAGIPnon}\\\t\${uniqSpeciesFLAGIP}\\\t\${uniqReadsFLAGIPnon}\\\t\${uniqReadsFLAGIP}\" >> ${OUTPUTDIR}/stat.log " >>${parafile}
		done
	done
done
[ $? == 0 ] && \
if [[ ! -f ${parafile}.completed ]] || [[ -f $parafile.failed_commands ]]
then
	CPUN=`wc -l $parafile |cut -f1 -d" "` && \
	ParaFly -c $parafile -CPU 1 -failed_cmds $parafile.failed_commands
fi
[ $? == 0 ] && \
touch ${OUTPUTDIR}/.status.${STEP}.IPuniqBound
STEP=$((STEP+1))