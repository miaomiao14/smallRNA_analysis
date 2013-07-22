#!/bin/bash -x


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
export PATH=${PIPELINE_DIRECTORY}/bin:$PATH




#########
# USAGE #
#########
# usage function
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -i input_file.[norm].bed -t type -o output_directory[current directory] -c cpu[8] 

This is a small RNA downstream analysis pipeline developed in the Zamore Lab in 
University of Massachusetts Medical School. 
Please email wei.wang2@umassmed.edu for any questions or bugs. 
Thanks for using it. 

OPTIONS:
	-h      Show this message
	-i      Input file in bed/norm.bed format, with full directory
	-o      Output directory, default: current directory
	-t      input file format:bed or normbed
	-c      Number of CPUs to use, default: 8
	-f      Configure file to provide/overwrite variables	

EOF
echo -en "\e[0m"
}
# taking options
while getopts "hi:c:o:t:f:" OPTION
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

FILE=${INPUT##*/}
piRNA_distance_distribution.pl $INPUT $TYPE $OUTDIR
Rscript piRNA_distance_plot.r ${OUTDIR}/${FILE}.5-5.distance.distribution
#Rscript piRNA_distance_plot.r 
