#! /bin/bash -x

################
# Major Config #
################
export PIPELINE_DIRECTORY=/home/hanb/nearline/small_RNA_Pipeline
export PATH=${PIPELINE_DIRECTORY}/bin:$PATH

#########
# USAGE #
#########
# usage function
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -u unox_folder -o ox_folder -O output_directory -s mouse|fly

This script takes the small RNA pipeline output folders of UNOX and OX library
and calculate the relative depth using edgeR.

Output will be in current folder.

OPTIONS:
	-h      Show this message
	-u      Directory of small RNA pipeline Unox library
	-o      Directory of small RNA pipeline OX library
	-s      Species [  mouse | fly ]

EOF
echo -en "\e[0m"
}

# taking options
while getopts "hu:o:s:O:" OPTION
do
	case $OPTION in
		h)
			usage && exit 1
		;;
		u)
			UNOX_DIR=`readlink -f $OPTARG`
		;;
		o)
			OX_DIR=`readlink -f $OPTARG`
		;;
		s)
			ORGANISM=$OPTARG
		;;
		O)
			OUTDIR=$OPTARG
		;;
		?)
			usage && exit 1
		;;
	esac
done

# if FQ or ORGANISM is undefined, print out usage and exit
if [[ -z $UNOX_DIR ]] || [[ -z $ORGANISM ]] || [[ -z $OX_DIR ]] 
then
	usage && exit 1
fi

# test existence
[ ! -d $UNOX_DIR ] && echo -e "\e[1;31mError: $UNOX_DIR does not exist...\e[0m" && exit 1
[ ! -d $OX_DIR ]   && echo -e "\e[1;31mError: $OX_DIR does not exist...\e[0m" && exit 1
PREFIX=`echo -e "${UNOX_DIR##*/}\n${OX_DIR##*/}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && PREFIX=${PREFIX%.*}
[ -z "${PREFIX}" ] && PREFIX=${UNOX_DIR##*/}"_OX"
[ ! -z $OUTDIR ] || OUTDIR=$PWD
# test wether we need to create the new directory
mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the direcory of current folder\e[0m"
# enter destination direcotry, since later it will be ran in this directory, there is NO need to add OUTPUT DIR in the following codes
cd ${OUTDIR} || echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" || exit 1
# test writtability
touch ${OUTDIR}/.writting_permission && rm -rf ${OUTDIR}/.writting_permission || echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" || exit 1

# check finish status
UNOX_VERSION=`ls -a $UNOX_DIR | grep .finish`
OX_VERSION=`ls -a $OX_DIR | grep .finish`
# check version of each run
[ $UNOX_VERSION != $OX_VERSION ] && \
	echo -e "\e[1;31mError:\n(1) ${UNOX_VERSION#*finished_}\n(2) ${OX_VERSION#*finished_}\nError: $UNOX_DIR and $OX_DIR were not ran by the same version of small RNA pipeline\e[0m" && \
	exit 1

case ${ORGANISM} in
mouse)
	GENOME=mm9
	piRNA_LEN=23
;;
fly)
	GENOME=dm3
	piRNA_LEN=22
;;
*)
	echo -e "\e[1;31mError: Unknown orgnanism... Currently only mouse/fly is supported\e[0m"
	exit 2
;;
esac

UNOX_ALIGNED_INSERT_FILE=${UNOX_DIR}/`ls $UNOX_DIR | grep $GENOME | grep al.insert`
[ ! -f $UNOX_ALIGNED_INSERT_FILE ] && \
 echo -e "\e[1;31mError: $UNOX_DIR does not have al.insert file, please make sure that you didn't not delete it...\e[0m" && exit 1

OX_ALIGNED_INSERT_FILE=${OX_DIR}/`ls $OX_DIR | grep $GENOME | grep al.insert`
[ ! -f $OX_ALIGNED_INSERT_FILE ] && \
 echo -e "\e[1;31mError: $OX_DIR does not have al.insert file, please make sure that you didn't not delete it...\e[0m" && exit 1

# merge
[ ! -f .status.${PREFIX}.mergeInsert ] && \
	myMerge -i $UNOX_ALIGNED_INSERT_FILE $OX_ALIGNED_INSERT_FILE -o ${PREFIX}.Unox_Ox.al.insert.merged -c 1 -t 2 && \
touch .status.${PREFIX}.mergeInsert

awk -v minlen=piRNA_LEN $'BEGIN{OFS="\t"}{if ($2>0 && $3>0 && length($1)>minlen) {print $0}}' ${PREFIX}.Unox_Ox.al.insert.merged > ${PREFIX}.Unox_Ox.al.insert.merged.noZero && \
Rscript $PIPELINE_DIRECTORY/bin/calculateUnoxOxEffectiveDepthFromMergedInsert.R ${PREFIX}.Unox_Ox.al.insert.merged.noZero > ${PREFIX}.Unox_Ox.al.insert.merged.noZero.depth && \
awk 'BEGIN{getline} {FS=" "; unoxSize=$3*$4;getline;oxSize=$3*$4; printf "%.3f", oxSize/unoxSize}' ${PREFIX}.Unox_Ox.al.insert.merged.noZero.depth >  ${PREFIX}.Unox_Ox.al.insert.merged.noZero.multiFactor
