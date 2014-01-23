#!/bin/bash


################
# Major Config #
################
# pipeline version
export smallRNA_downstream_analysis=1.0.0
# this pipeline is still in debug mode
export DEBUG=1 
# pipeline address: if you copy all the files to another directory, this is the place to change; under this directory sits two directories, bin and common. bin stores all the binary executables and common stores all the information of each ORGANISM.
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
# set PATH to be aware of pipeline/bin; this bin directory needs to be searched first
export PATH=${PIPELINE_DIRECTORY}/:$PATH




#########
# USAGE #
#########
# usage function
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -i input_file.[norm].bed -t type -o output_directory[current directory] 

This is a small RNA downstream analysis pipeline. 

OPTIONS:
	-h      Show this message
	-i      Input file in bed/norm.bed format, with full directory
	-o      Output directory, default: current directory
	-t      input file format:bed or normbed
	-c      Number of CPUs to use, default: 8
	-r		SRA or DEG
	-f      Configure file to provide/overwrite variables	

EOF
echo -en "\e[0m"
}
# taking options
while getopts "hi:c:o:t:f:r:" OPTION
do
	case $OPTION in
		h)
			usage && exit 1
		;;
		i)
			INPUT=$OPTARG
		;;
		o)
			OUTDIR=$OPTARG
		;;
		c)
			CPU=$OPTARG
		;;
		t)
			TYPE=$OPTARG
		;;
		r)
			RNA=$OPTARG
		;;
		f)
			CONFIG_FILE=$OPTARG
		;;
		?)
			usage && exit 1
		;;
	esac
done
# if FQ or ORGANISM is undefined, print out usage and exit
if [[ -z $INPUT ]] || [[ -z $TYPE ]] 
then
	usage && exit 1
fi
[ ! -z $OUTDIR ] || OUTDIR=$PWD
FILE=${INPUT##*/}
if [ $RNA = "DEG" ]
then
${PIPELINE_DIRECTORY}/piRNA_distance_distribution_DEG.pl $INPUT $TYPE $OUTDIR
else
${PIPELINE_DIRECTORY}/piRNA_distance_distribution_SRA.pl $INPUT $TYPE $OUTDIR
fi

${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/piRNA_distance_plot.r plot_distribution_summary ${OUTDIR}/${FILE}.5-5.distance.distribution.summary
${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/piRNA_distance_plot.r plot_distribution ${OUTDIR}/${FILE}.5-5.distance.distribution
