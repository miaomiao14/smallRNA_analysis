#!/bin/bash -x
# small RNA pipeline in the Zamore Lab
# single library mode
# 2014-1-14
# Bo W Han (bo.han@umassmed.edu)
# Phillip Zamore Lab
# RNA Therapeutics Institute
# HHMI & University of Massachusetts Medical School

################
# Major Config #
################
# pipeline version
export smallRNA_Pipeline_Version=1.6.4
# this pipeline is still in debug mode
export DEBUG=1 
# pipeline address: if you copy all the files to another directory, this is the place to change; under this directory sits two directories, bin and common. bin stores all the binary executables and common stores all the information of each ORGANISM.
export PIPELINE_DIRECTORY=/home/hanb/nearline/small_RNA_Pipeline
# set PATH to be aware of pipeline/bin; this bin directory needs to be searched first
export PATH=${PIPELINE_DIRECTORY}/bin:$PATH
# for python
export PYTHONPATH=/home/hanb/lib/python2.7/site-packages:${PYTHONPATH}
# for C++ run time linking
export LD_LIBRARY_PATH=/home/hanb/lib64:/home/hanb/lib:${LD_LIBRARY_PATH}

#########
# USAGE #
#########
# usage function
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -i input_file.fq[.gz] -s fly|mouse -o output_directory[current directory] -c cpu[8] 

Please pay attention to the version number.

This is a single library small RNA pipeline developed in the Zamore Lab in 
University of Massachusetts Medical School. 

Please email Bo.Han@umassmed.edu for any questions or bugs. 

Thanks for using it. 

OPTIONS:
	-h      Show this message
	-i      Input file in fastq or gzip -fped fastq format, with full directory
	-o      Output directory, default: current directory
	-s      Species [ human | mouse | mm10 | fly ]
	-c      Number of CPUs to use, default: 8
	-f      Configure file to provide/overwrite variables	

EOF
echo -en "\e[0m"
}
# taking options
while getopts "hi:c:o:s:f:" OPTION
do
	case $OPTION in
		h)
			usage && exit 1
		;;
		i)
			FQ=$OPTARG # FQ now is the string provided by the user
		;;
		o)
			OUTDIR=`readlink -f $OPTARG` # OUTDIR now is the absolute path
		;;
		c)
			CPU=$OPTARG
		;;
		s)
			ORGANISM=$OPTARG
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
if [[ -z $FQ ]] || [[ -z $ORGANISM ]] 
then
	usage && exit 1
fi

###########################
# running necessity check #
###########################
# function to check whether current shell can find all the software/programs needed to finish running the pipeline
function checkExist {
	echo -ne "\e[1;32m\"${1}\" is using: \e[0m" && which "$1"
	[[ $? != 0 ]] && echo -e "\e[1;36mError: cannot find software/function ${1}! Please make sure that you have installed the pipeline correctly.\nExiting...\e[0m" && \
	exit 1
}
echo -e "\e[1;35mTesting required softwares/scripts:\e[0m"
checkExist "echo" # used to print message
checkExist "date" # used to get time
checkExist "tr" # used to convert between up/lower cases
checkExist "mkdir" # used to creat new directory
checkExist "readlink" # used to find the real location of a file
checkExist "rm" # used to delete file/directory
checkExist "mv" # used to move or rename file
checkExist "sort" # used to sort alphnum
checkExist "touch" # used to create new empty file
checkExist "md5sum" # used to create a uid so that the status file for different run will not mess up eath other
checkExist "seq" # used to create serial number
checkExist "awk" # a programming language
checkExist "grep" # a programming language
checkExist "python" # a programming language
checkExist "ruby" # a programming language
checkExist "samtools" # tools to process sam/bam file
checkExist "bedtools" # tools to process bed formatted file
checkExist "gs" # used to merge pdf
checkExist "gzip" # used to compress file
checkExist "Rscript" # a programming language
checkExist "bowtie" # mapping tool, used for most mapping purpose
checkExist "bowtie-build" # building index for bowtie
checkExist "gsnap" # mapping tool, used for mapping to exon-exon junction
checkExist "gtf_splicesites" # used to create splice site file for gsnap
checkExist "iit_store" # used to create splice site file for gsnap
checkExist "genomeCoverageBed2" # modified bedtools to take in-house bed2 format and make bedgraph format; !! it doesn't consider the 5th column and only take integer as the 4th col, thus should only be used for unnormalized unique mappers
checkExist "bedItemOverlapCount2" # modified bedItemOverlapCount from kent source to take in-house bed2 format and make bedgrah format; it can handle NTM but doesn't have the option -strand and -scale
checkExist "weblogo" # make weblogo 
checkExist "information_content.py" # calculate information content from fasta
checkExist "ParaFly" # a handy tool in trinity package to wrap a serials independent commands and run them with multiple threads
checkExist "bedGraphToBigWig" # a tool from kent source to convert bedgraph to bigWig format
checkExist "phred_test" # in-house C++ program to find out wether a fastq format is in phred+33 or phred+64 mode
checkExist "ppbed2" # in-house C++ program to calculate ping-pong score from two bed2 files
checkExist "ppbed2each" # in-house C++ program to calculate ping-pong score of each chr/transcript from two bed2 file
checkExist "istBedToBed2" # in-house C++ program to convert bed file, which is made by bowtie-mapping of insert file to bed2 format
checkExist "calculate_miRNA_heterogeneity" # in-house C++ program to calculate miRNA heterogeneity from bed2 format
checkExist "fastq2insert" # in-house C++ program to convert compressed/uncompressed fastq file to insert file
checkExist "direct_mapping.sh" # direct mapping script
echo -e "\e[1;35mDone with testing required softwares/scripts, starting pipeline...\e[0m"

##################
# File/Dir check #
##################
# check whether input file exist
[ ! -f $FQ ] && echo -e "\e[1;31mError: Cannot find input file $FQ\e[0m" && exit 1
# get full path of the FQ and use the full path when getting the insert file.
FULLPATH_FQ=`readlink -f $FQ`
# use $FQ later only as filename
FQ=`basename $FULLPATH_FQ`
PREFIX=${FQ%.f[qa]*}
# if CPU is undefined or containing non-numeric char, then use 8
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8

# OUTPUT directory defined? otherwise use current directory
[ ! -z $OUTDIR ] || OUTDIR=$PWD
# test wether we need to create the new directory
mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the direcory of input fastq file\e[0m"
# enter destination direcotry, since later it will be ran in this directory, there is NO need to add OUTPUT DIR in the following codes
cd ${OUTDIR} || (echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" && exit 1)
# test writtability
touch .writting_permission && rm -rf .writting_permission || (echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" && exit 1)
# directory storing the pdfs
PDF_DIR=pdfs && mkdir -p $PDF_DIR

#############
# Variables #
#############
# step counter
STEP=1
# unique job id
JOBUID=`echo ${FQ##*/} | md5sum | cut -d" " -f1`
# formating date for log
ISO_8601='%Y-%m-%d %H:%M:%S %Z'
# table to store the basic statistics of the library (genomic mappability). PS: counts of each genomic features will be stored in .summary file, produced by intersect_all.sh 
TABLE=${PREFIX}.basic_stats
# log file to store the timing
LOG=${PREFIX}.pipeline_${smallRNA_Pipeline_Version}.log
# set ORGANISM, convert upper cases into lower case so it is easier to write the case statement later
ORGANISM=`echo ${ORGANISM} | tr '[A-Z]' '[a-z]'`
# assign different values to the generalized variables (same name for different ORGANISMs) according to which ORGANISM fed
case ${ORGANISM} in
human) #TODO
	echo "not yet" && exit 1
	Genome=hg19
	COMMON_FOLDER=$PIPELINE_DIRECTORY/common_files/$Genome
	rRNA_MM=2
	hairpin_MM=1
	genome_MM=1
	cluster_MM=1
	transposon_MM=3
	mRNA_MM=1
;;
mouse)
# using mm9 for mouse 
	Genome=mm9
	COMMON_FOLDER=$PIPELINE_DIRECTORY/common_files/$Genome
	rRNA_MM=3
	hairpin_MM=1
	genome_MM=1
	cluster_MM=1
	transposon_MM=3
	mRNA_MM=1
;;
fly)
# using dm3 for fly
	Genome=dm3
	COMMON_FOLDER=$PIPELINE_DIRECTORY/common_files/$Genome
	rRNA_MM=3
	hairpin_MM=0
	genome_MM=0
	cluster_MM=1
	transposon_MM=2
	mRNA_MM=0
;;
*)
	echo -e "\e[1;31mError: Unknown orgnanism... Currently only mouse/fly is supported\e[0m"
	exit 2
;;
esac

# fasta file for the genome
GENOME_FA=$COMMON_FOLDER/genome.fa
# bowtie1 index for the genome
BOWTIEINDEX_Genome=$COMMON_FOLDER/$Genome
# bowtie1 index for rRNA, from iGenome
BOWTIEINDEX_rRNA=$COMMON_FOLDER/ribosomal
# bowtie1 index for the hairpin (miRNA) of this ORGANISM, from miRBase
BOWTIEINDEX_Hairpin=$COMMON_FOLDER/hairpin
# bowtie1 index for the hairpin (miRNA) of this ORGANISM, from miRBase
BOWTIEINDEX_Mature=$COMMON_FOLDER/mature
# chrom information of this ORGANISM
CHROM=$COMMON_FOLDER/ChromInfo.txt
# gtf file (required by gsnap use map intron)
GTF=$COMMON_FOLDER/genes.gtf

#######>>> Begin of custumer variables <<<#######
# PreMappingList stores an array of variable names. those variables names store the address of bowtie indexes. 
# this pipeline will map fastq to each one of those bowite indexes sequentially
# mapped reads will be excluded from mapping to next index
# premappinglist mapping happens after rRNA and miRNA hairpin mapping
case $ORGANISM in 
human)
	declare -a PreMappingList=(" ") 
	declare -a PreMappingMM=(" ")
;;
mouse)
	mmu_virus=$COMMON_FOLDER/mmu_virus
	declare -a PreMappingList=("mmu_virus") 
	declare -a PreMappingMM=("2")
#----------start_of_example----------
# how to add new targets
	#new_target=$COMMON_FOLDER/mmu_virus
	#PreMappingList=("${PreMappingList[@]}" "new_target")
	#PreMappingMM=("${PreMappingMM[@]}" "2")
#----------end_of_example----------
;;
fly)
	dmel_virus=$COMMON_FOLDER/dmel_virus
	declare -a PreMappingList=("dmel_virus") 
	declare -a PreMappingMM=("2")	
#----------start_example----------
# how to add new targets
	#new_target=$COMMON_FOLDER/dmel_virus
	#PreMappingList=("${PreMappingList[@]}" "new_target")
	#PreMappingMM=("${PreMappingMM[@]}" "2")
#----------end_of_example----------
;;
esac
#######>>> End of custumer variables <<<########

#######################
# Function definition #
#######################
# count reads for bed2 format
function bedwc {
	awk '{a[$7]=$4}END{COUNTER=0; for (b in a){COUNTER+=a[b]} printf "%d" , COUNTER}' $1
}
export -f bedwc 
# produce length distribution for bed2 format
function bed2lendis {
	awk '{l=$3-$2; if (l>m) m=l; c[l]+=$4/$5;}END{for (d=1;d<=m;++d) {printf "%d\t%.0f\n", d, (c[d]?c[d]:0)} }' $1
}
export -f bed2lendis

# reading config file
[ ! -z $CONFIG_FILE ] && [ -f $CONFIG_FILE ] && \
. $CONFIG_FILE && \
echo -e "Reading in configure for custom variable assignment..." || \
echo -e "Warning: Failed to source configure file $CONFIG_FILE"

##############################
# beginning running pipeline #
##############################
# if log file already exists, then "resume" it. otherwise "begin" it.
[ ! -f $LOG ] && \
echo -e "`date "+$ISO_8601"`\tbeginning running Zamore Lab small RNA pipeline version $smallRNA_Pipeline_Version in one lib mode"   | tee -a $LOG || \
echo -e "`date "+$ISO_8601"`\t...resuming running Zamore Lab small RNA pipeline version $smallRNA_Pipeline_Version in one lib mode" | tee -a $LOG 

########################################
## Pre Processing before any Mapping ###
########################################
# convering fastq to insert; quality information will be lost
echo -e "`date "+$ISO_8601"`\tconverting fastq format into insert format" | tee -a $LOG
READS_DIR=input_read_files && mkdir -p $READS_DIR
# insert file, a format with two fields delimited by a tab. Sequence and number of times it was read, used to save time/space; quality information is lost
INSERT=$READS_DIR/${PREFIX}.insert
[ ! -f .${JOBUID}.status.${STEP}.fq2insert ] && \
	fastq2insert ${FULLPATH_FQ} ${INSERT} && \
	touch .${JOBUID}.status.${STEP}.fq2insert 
STEP=$((STEP+1))
	
#####################################
# Pre Processing before any Mapping #
#####################################
# getting rid of sequences mapping to rRNA with 2 mismatches (by default) allowed
## as soon as one read is mappable under -v 2 condition, it is dumped, so -k 1 is used for speed purpose, and no need to be "best" 
echo -e "`date "+$ISO_8601"`\tmapping to rRNA, with $rRNA_MM mismatch(es) allowed" | tee -a $LOG
rRNA_DIR=rRNA_mapping && mkdir -p $rRNA_DIR
rRNA_BED2=$rRNA_DIR/${PREFIX}.rRNA.bed2
x_rRNA_INSERT=$READS_DIR/${PREFIX}.x_rRNA.insert
[ ! -f .${JOBUID}.status.${STEP}.rRNA_mapping ] && \
	bowtie -r -S -v $rRNA_MM -k 1 -p $CPU \
		--un $x_rRNA_INSERT \
		$BOWTIEINDEX_rRNA \
		${INSERT} \
		1> /dev/stdout \
		2> /dev/null | \
	samtools view -bSF 0x4 - 2>/dev/null | \
	bedtools bamtobed -i - > ${PREFIX}.rRNA.insertBed && \
	istBedToBed2 ${INSERT} ${PREFIX}.rRNA.insertBed > $rRNA_BED2 && \
	rm -rf ${PREFIX}.rRNA.insertBed && \
    touch .${JOBUID}.status.${STEP}.rRNA_mapping # by checking all the status, ll -a | grep status, you will know which step has finished 
totalReads=`awk '{a+=$2}END{printf "%d", a}' ${INSERT}`
rRNAReads=`bedwc $rRNA_BED2`
nonrRNAReads=$((totalReads-rRNAReads))
## increment the step counter, no matter previous step succeed or not. Using a variable will make it easy to insert/delete steps
STEP=$((STEP+1))
    
#########################
# miRNA hairpin Mapping #
#########################
# mapping to miRNA hairpins
echo -e "`date "+$ISO_8601"`\tmapping to hairpin, with $hairpin_MM mismatch(es) allowed" | tee -a $LOG
# dir to store hairpin-based mapping
MIRNA_DIR=hairpins_mapping && mkdir -p $MIRNA_DIR
# dir to store genome-based mapping
GENOMIC_MAPPING_DIR=genome_mapping && mkdir -p $GENOMIC_MAPPING_DIR
# hairpin mappers in insert format
x_rRNA_HAIRPIN_INSERT=$READS_DIR/${PREFIX}.x_rRNA.hairpin.insert
# non-hairpin mapper in insert format
x_rRNA_x_hairpin_INSERT=$READS_DIR/${PREFIX}.x_rRNA.x_hairpin.insert
# hairpin mapping in bed2 format
x_rRNA_HAIRPIN_BED2=$MIRNA_DIR/${PREFIX}.x_rRNA.hairpin.v${hairpin_MM}a.bed2
x_rRNA_HAIRPIN_BED2_LENDIS=$MIRNA_DIR/${PREFIX}.x_rRNA.hairpin.v${hairpin_MM}a.lendis
# genome mapping of hairpin mappers in bed2 format
x_rRNA_HAIRPIN_GENOME_BED2=$GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed2
x_rRNA_HAIRPIN_GENOME_LOG=$GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA.hairpin.${Genome}v${genome_MM}a.log

[ ! -f .${JOBUID}.status.${STEP}.hairpin_mapping ] && \
	bowtie -r -v $hairpin_MM -a --best --strata -p $CPU -S \
		${BOWTIE_PHRED_OPTION} \
		--al $x_rRNA_HAIRPIN_INSERT \
		--un $x_rRNA_x_hairpin_INSERT \
		$BOWTIEINDEX_Hairpin \
		$x_rRNA_INSERT \
		1> /dev/stdout \
		2> /dev/null  | \
	samtools view -bSF 0x4 - 2>/dev/null | \
	bedtools bamtobed -i - | awk '$6=="+"' > ${PREFIX}.x_rRNA.hairpin.v${hairpin_MM}a.bed && \
	istBedToBed2 $x_rRNA_INSERT ${PREFIX}.x_rRNA.hairpin.v${hairpin_MM}a.bed > $x_rRNA_HAIRPIN_BED2 && \
	rm -rf ${PREFIX}.x_rRNA.hairpin.v${hairpin_MM}a.bed && \
	bed2lendis $x_rRNA_HAIRPIN_BED2 > $x_rRNA_HAIRPIN_BED2_LENDIS && \
	bowtie -r -v $genome_MM -a --best --strata -p $CPU \
		-S \
		$BOWTIEINDEX_Genome \
		$x_rRNA_HAIRPIN_INSERT \
		1> /dev/stdout \
		2> $x_rRNA_HAIRPIN_GENOME_LOG | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	bedtools bamtobed -i - > ${PREFIX}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed && \
	istBedToBed2 $x_rRNA_HAIRPIN_INSERT ${PREFIX}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed > $x_rRNA_HAIRPIN_GENOME_BED2 && \
	rm -rf ${PREFIX}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed && \
	touch .${JOBUID}.status.${STEP}.hairpin_mapping
STEP=$((STEP+1))
hairpinReads=`bedwc $x_rRNA_HAIRPIN_BED2`
# run miRNA pipeline
echo -e "`date "+$ISO_8601"`\trunning miRNA analysis pipeline" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.miRNA_pipeline ] && \
	miRNA_pipeline.sh $BOWTIEINDEX_Hairpin $BOWTIEINDEX_Mature $x_rRNA_HAIRPIN_BED2  && \
	touch .${JOBUID}.status.${STEP}.miRNA_pipeline
STEP=$((STEP+1))

######################
# Pre Genome Mapping #
######################
# count how many indexed need to map
COUNTER=0
# running variable
INPUT=$x_rRNA_x_hairpin_INSERT
# if there are any indexes need to be run
[[ ${#PreMappingList[@]} > 0 ]] && \
for COUNTER in `seq 0 $((${#PreMappingList[@]}-1))`; do \
	TARGET=${PreMappingList[$COUNTER]} && \
	OUTDIR1=${TARGET}_mapping && mkdir -p $OUTDIR1 && \
	PREFIX1=`basename $INPUT` && PREFIX1=${PREFIX1%.insert} && \
	MM=${PreMappingMM[$COUNTER]} && \
	echo -e "`date "+$ISO_8601"`\tmapping to ${TARGET}, with $MM mismatch(es) allowed" | tee -a $LOG && \
	[ ! -f .${JOBUID}.status.${STEP}.${TARGET}_mapping ] && \
		bowtie -r -v $MM -a --best --strata -p $CPU -S \
			${BOWTIE_PHRED_OPTION} \
			--un ${INPUT%.insert}.x_${TARGET}.insert \
			${!TARGET} \
			$INPUT \
			1> /dev/stdout \
			2> /dev/null | \
		samtools view -bSF 0x4 - 2>/dev/null | bedtools bamtobed -i - > ${PREFIX1}.${TARGET}.v${hairpin_MM}a.bed && \
		istBedToBed2 $INPUT ${PREFIX1}.${TARGET}.v${hairpin_MM}a.bed > $OUTDIR1/${PREFIX1}.${TARGET}.v${hairpin_MM}a.bed2 && \
		rm -rf ${PREFIX1}.${TARGET}.v${hairpin_MM}a.bed && \
		touch .${JOBUID}.status.${STEP}.${TARGET}_mapping
	# mappingReads=`bedwc ${INPUT%.f[qa]*}.${TARGET}.v${hairpin_MM}a.bed2` # currently not used
	STEP=$((STEP+1))
	# set OUPUT of this loop to the INPUT of next
	INPUT=${INPUT%.insert}.x_${TARGET}.insert
done

##################
# Genome Mapping #
##################
# take the OUTPUT of last step as INPUT
INSERT=`basename ${INPUT}`
# bed2 format storing all mappers for genomic mapping
GENOME_ALLMAP_BED2=$GENOMIC_MAPPING_DIR/${INSERT%.insert}.${Genome}v${genome_MM}.all.bed2
GENOME_ALLMAP_LOG=$GENOMIC_MAPPING_DIR/${INSERT%.insert}.${Genome}v${genome_MM}.all.log

# bed2 format storing all mappers by gsnap
GSNAP_ALLMAP_BED12=$GENOMIC_MAPPING_DIR/${INSERT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed12
GSNAP_ALLMAP_BED6=$GENOMIC_MAPPING_DIR/${INSERT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed6
GSNAP_ALLMAP_LOG=$GENOMIC_MAPPING_DIR/${INSERT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.log

# bed2 format storing unique mappers for genomic mapping
GENOME_UNIQUEMAP_BED2=$GENOMIC_MAPPING_DIR/${INSERT%.insert}.${Genome}v${genome_MM}.unique.bed2

# bed2 format storing unique mappers for genomic mapping and miRNA hairpin mapper
GENOME_UNIQUEMAP_BED2_hairpin=$GENOMIC_MAPPING_DIR/${INSERT%.insert}.${Genome}v${genome_MM}.unique.+hairpin.bed2
# mapping insert file to genome
echo -e "`date "+$ISO_8601"`\tmapping to genome, with ${genome_MM} mismatch(es) allowed " | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.genome_mapping ] && \
	bowtie -r -v $genome_MM -a --best --strata -p $CPU \
		--al  ${INPUT%.insert}.${Genome}v${genome_MM}a.al.insert \
		--un  ${INPUT%.insert}.${Genome}v${genome_MM}a.un.insert \
		-S \
		$BOWTIEINDEX_Genome \
		${INPUT} \
		1> /dev/stdout \
		2> $GENOME_ALLMAP_LOG | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	bedtools bamtobed -i - > ${INSERT%.insert}.${Genome}v${genome_MM}a.insert.bed && \
	istBedToBed2 $INPUT ${INSERT%.insert}.${Genome}v${genome_MM}a.insert.bed > ${GENOME_ALLMAP_BED2} && \
	rm -rf ${INSERT%.insert}.${Genome}v${genome_MM}a.insert.bed && \
	touch .${JOBUID}.status.${STEP}.genome_mapping
STEP=$((STEP+1))

# produce iit format for intron annotation used in gsnap, if it doesn't exist
[ ! -f ${COMMON_FOLDER}/${Genome}/${Genome}.splicesites.iit ] && \
	cat $GTF | gtf_splicesites > ${COMMON_FOLDER}/${Genome}/${Genome}.splicesites && \
	cat ${COMMON_FOLDER}/${Genome}/${Genome}.splicesites | iit_store -o ${COMMON_FOLDER}/${Genome}/${Genome}.splicesites.iit
# converting non-mapped insert file to fasta for gsnap mapping
echo -e "`date "+$ISO_8601"`\tmapping unmapped read to genome using gsnap [warning: gsnap cannot handle reads less than 17nt]" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.gsnap_genome_mapping ] && \
	awk '{print ">"$1"_"$2"\n"$1}' ${INPUT%.insert}.${Genome}v${genome_MM}a.un.insert > ${INPUT%.insert}.${Genome}v${genome_MM}a.un.fa && \
	gsnap -A sam -D $COMMON_FOLDER -d $Genome \
		-m ${genome_MM} \
		--maxsearch=1000000 \
		-n 1000000 \
		--quiet-if-excessive \
		-t $CPU \
		-N 0 \
		-i 100 \
		-s ${COMMON_FOLDER}/${Genome}/${Genome}.splicesites.iit \
		${INPUT%.insert}.${Genome}v${genome_MM}a.un.fa \
		2> $GSNAP_ALLMAP_LOG | \
	samtools view -bS -F0x4 - 2>/dev/null | bedtools bamtobed -bed12 -i - > ${GSNAP_ALLMAP_BED12}.temp && \
	awk 'BEGIN{OFS="\t"}NR==FNR{++ct[$4]} {$5=ct[$4]; split($4,arr,"_"); $4=arr[2];print}' ${GSNAP_ALLMAP_BED12}.temp ${GSNAP_ALLMAP_BED12}.temp > ${GSNAP_ALLMAP_BED12} && \
	awk '{split($4,arr,"_"); ct[$4]=arr[2]}END{for (a in ct){total+=ct[a]} print total}' ${GSNAP_ALLMAP_BED12}.temp > ${GSNAP_ALLMAP_BED12}.count && \
	awk '$5==1&&$10>1' ${GSNAP_ALLMAP_BED12} | bedtools bed12tobed6 -i stdin > ${GSNAP_ALLMAP_BED6} && \
	rm -rf ${INPUT%.insert}.${Genome}v${genome_MM}a.un.fa ${GSNAP_ALLMAP_BED12}.temp ${GSNAP_ALLMAP_BED12} && \
	touch .${JOBUID}.status.${STEP}.gsnap_genome_mapping
spliceMapTotal=`cat ${GSNAP_ALLMAP_BED12}.count`
STEP=$((STEP+1))

# separating unique and multiple mappers
echo -e "`date "+$ISO_8601"`\tseparating unique and multiple mappers" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.separate_unique_and_multiple ] && \
	awk 'BEGIN{OFS="\t"}{if ($5==1) print $0}' ${GENOME_ALLMAP_BED2} \
	1> ${GENOME_UNIQUEMAP_BED2}	&& \
	touch .${JOBUID}.status.${STEP}.separate_unique_and_multiple
totalMapCount=`bedwc ${GENOME_ALLMAP_BED2}`
uniqueMapCount=`bedwc ${GENOME_UNIQUEMAP_BED2}`
multipMapCount=$((totalMapCount-uniqueMapCount))
STEP=$((STEP+1))

#####################
# Length Separation #
#####################
# length separation
case ${ORGANISM} in
mouse)
	echo -e "`date "+$ISO_8601"`\tseparating reads based on length" | tee -a $LOG
	[ ! -f .${JOBUID}.status.${STEP}.sep_length ] && \
		paraFile=${RANDOM}${RANDOM}.para && \
		echo "awk '\$3-\$2==19' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_BED2%bed2}19mer.bed2" >  $paraFile && \
		echo "awk '\$3-\$2>=20 && \$3-\$2<23' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_BED2%bed2}20-22mer.bed2" >> $paraFile && \
		echo "awk '\$3-\$2>=23' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_BED2%bed2}L23mer.bed2" >> $paraFile && \
		ParaFly -c $paraFile -CPU $CPU && \
		rm -rf ${paraFile}* && \
		touch .${JOBUID}.status.${STEP}.sep_length
	  STEP=$((STEP+1))
;;
## fly specific analysis
fly)
	echo -e "`date "+$ISO_8601"`\tseparating reads based on length" | tee -a $LOG
	[ ! -f .${JOBUID}.status.${STEP}.sep_length ] && \
		paraFile=${RANDOM}${RANDOM}.para && \
		echo "awk '\$3-\$2< 20' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_BED2%bed2}l20mer.bed2" > $paraFile && \
		echo "awk '\$3-\$2>=20 && \$3-\$2<23' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_BED2%bed2}20-22mer.bed2" >> $paraFile && \
		echo "awk '\$3-\$2>=23' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_BED2%bed2}L23mer.bed2" >> $paraFile && \
		ParaFly -c $paraFile -CPU $CPU && \
		rm -rf ${paraFile}* && \
		touch  .${JOBUID}.status.${STEP}.sep_length
	STEP=$((STEP+1))
;;
*)
	echo "[$LINENO] Error: unrecognized organism $ORGANISM"
	exit 2
;;
esac

# plotting length distribution
echo -e "`date "+$ISO_8601"`\tmaking length distribution" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.plotting_length_dis ] && \
	awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}' ${GENOME_ALLMAP_BED2}  | sort -k1,1n > ${GENOME_ALLMAP_BED2}.lendis && \
	awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}' ${GENOME_UNIQUEMAP_BED2}  | sort -k1,1n > ${GENOME_UNIQUEMAP_BED2}.lendis && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_lendis.R ${GENOME_ALLMAP_BED2}.lendis $PDF_DIR/`basename ${GENOME_ALLMAP_BED2%.bed*}`.x_hairpin && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_lendis.R ${GENOME_UNIQUEMAP_BED2}.lendis $PDF_DIR/`basename ${GENOME_UNIQUEMAP_BED2%.bed*}`.x_hairpin && \
	awk '{ct[$1]+=$2}END{for (l in ct) {print l"\t"ct[l]}}' ${GENOME_ALLMAP_BED2}.lendis $x_rRNA_HAIRPIN_BED2_LENDIS | sort -k1,1n > ${GENOME_ALLMAP_BED2}.+hairpin.lendis && \
	awk '{ct[$1]+=$2}END{for (l in ct) {print l"\t"ct[l]}}' ${GENOME_UNIQUEMAP_BED2}.lendis $x_rRNA_HAIRPIN_BED2_LENDIS | sort -k1,1n > ${GENOME_UNIQUEMAP_BED2}.+hairpin.lendis && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_lendis.R ${GENOME_ALLMAP_BED2}.+hairpin.lendis $PDF_DIR/`basename ${GENOME_ALLMAP_BED2%.bed*}`.+hairpin && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_lendis.R ${GENOME_UNIQUEMAP_BED2}.+hairpin.lendis $PDF_DIR/`basename ${GENOME_UNIQUEMAP_BED2%.bed*}`.+hairpin && \
	touch .${JOBUID}.status.${STEP}.plotting_length_dis
STEP=$((STEP+1))


##################
# Print to table #
##################

echo -e "total reads as input of the pipeline\t${totalReads}" > $TABLE && \
echo -e "rRNA reads with ${rRNA_MM} mismatches\t${rRNAReads}" >> $TABLE && \
echo -e "genome mapping reads (-rRNA; +miRNA_hairpin)\t$((totalMapCount+hairpinReads+spliceMapTotal))" >> $TABLE && \
echo -e "miRNA hairpin reads\t${hairpinReads}" >> $TABLE && \
echo -e "genome mapping reads (-rRNA; -miRNA_hairpin)\t$((totalMapCount+spliceMapTotal))" >> $TABLE && \
echo -e "genome mapping reads (-rRNA; -miRNA_hairpin; -exon_exon_junction_mapper)\t${totalMapCount}" >> $TABLE && \
echo -e "genome unique mapping reads (-rRNA; -miRNA_hairpi; -exon_exon_junction_mapper)\t${uniqueMapCount}" >> $TABLE && \
echo -e "genome multiple mapping reads (-rRNA; -miRNA_hairpi; -exon_exon_junction_mapper)\t${multipMapCount}" >> $TABLE && \
echo -e "exon-exon junction reads (unmappable to genome)\t${spliceMapTotal}" >> $TABLE && \

###############################
# Intersecting Genome Feature #
###############################

# running intersecting bed
echo -e "`date "+$ISO_8601"`\trunning intersectBed for all features" | tee -a $LOG
SUMMARY_DIR=summaries && mkdir -p $SUMMARY_DIR
[ ! -f .${JOBUID}.status.${STEP}.intersect_all ] && \
INTERSECT_DIR1=intersect_all_length && \
mkdir -p $INTERSECT_DIR1 && \
intersect_all.sh \
	${GENOME_ALLMAP_BED2} \
	$SUMMARY_DIR/`basename ${GENOME_ALLMAP_BED2%.bed2}.all_len.summary` \
	$ORGANISM \
	$GENOME_FA \
	$CPU \
	$INTERSECT_DIR1 && \
INTERSECT_DIR2=intersect_piRNA_length && \
mkdir -p $INTERSECT_DIR2 && \
intersect_all.sh \
	${GENOME_ALLMAP_BED2%bed2}L23mer.bed2 \
	$SUMMARY_DIR/`basename ${GENOME_ALLMAP_BED2%.bed2}.piRNA.summary` \
	$ORGANISM \
	$GENOME_FA \
	$CPU \
	$INTERSECT_DIR2 && \
INTERSECT_DIR3=intersect_siRNA_length && \
mkdir -p $INTERSECT_DIR3 && \
intersect_all.sh \
	${GENOME_ALLMAP_BED2%bed2}20-22mer.bed2 \
	$SUMMARY_DIR/`basename ${GENOME_ALLMAP_BED2%.bed2}.siRNA.summary` \
	$ORGANISM \
	$GENOME_FA \
	$CPU \
	$INTERSECT_DIR3 && \
	touch .${JOBUID}.status.${STEP}.intersect_all
STEP=$((STEP+1))


########################
# Makeing BigWig Files #
########################
# normalization factor, currently using unique genome mappers
NormScale=`echo $((uniqueMapCount)) | awk '{printf "%f",1000000.0/$1}'`
# make BW files
echo -e "`date "+$ISO_8601"`\tmaking bigWig file for genome browser" | tee -a $LOG
BW_OUTDIR=bigWig && mkdir -p $BW_OUTDIR
[ ! -f .${JOBUID}.status.${STEP}.make_bigWig ] && \
	cat $x_rRNA_HAIRPIN_GENOME_BED2 ${GENOME_UNIQUEMAP_BED2} > $GENOME_UNIQUEMAP_BED2_hairpin && \
	bed22bw.sh \
		${GENOME_UNIQUEMAP_BED2_hairpin} \
		${GSNAP_ALLMAP_BED6} \
		${CHROM} \
		${NormScale} \
		$CPU \
		$ORGANISM \
		$BW_OUTDIR && \
	rm -rf $GENOME_UNIQUEMAP_BED2_hairpin && \
	touch .${JOBUID}.status.${STEP}.make_bigWig
STEP=$((STEP+1))

#################################
# Mapping to features directly #
#################################

# run bowtie mapping to repBase directly
BOWTIE_INDEX=${COMMON_FOLDER}/repBase/transposon
[ ! -f .${JOBUID}.status.${STEP}.repBase_mapping ] && \
	mkdir -p repBaseMapping && \
	direct_mapping.sh \
		${BOWTIE_INDEX} \
		$READS_DIR/${INSERT} \
		${transposon_MM} \
		repBaseMapping \
		"repBase" \
		$CPU \
		$NormScale && \
	touch .${JOBUID}.status.${STEP}.repBase_mapping
STEP=$((STEP+1))

# run bowtie mapping to Jia/Jie transposon directly
BOWTIE_INDEX=${COMMON_FOLDER}/transposon/transposon
[ ! -f .${JOBUID}.status.${STEP}.transposon_mapping ] && \
	mkdir -p transposonMapping && \
	direct_mapping.sh \
		${BOWTIE_INDEX} \
		$READS_DIR/${INSERT} \
		${transposon_MM} \
		transposonMapping \
		"transposon" \
		$CPU \
		$NormScale && \
	touch .${JOBUID}.status.${STEP}.transposon_mapping
STEP=$((STEP+1))
						
# run bowtie mapping to piRNA cluster
BOWTIE_INDEX=${COMMON_FOLDER}/piRNAcluster/piRNAcluster
[ ! -f .${JOBUID}.status.${STEP}.cluster_mapping ] && \
	mkdir -p piRNAclusterMapping && \
	direct_mapping2.sh \
		${BOWTIE_INDEX} \
		$READS_DIR/${INSERT} \
		${cluster_MM} \
		piRNAclusterMapping \
		"piRNAcluster" \
		$CPU \
		$NormScale && \
	touch .${JOBUID}.status.${STEP}.cluster_mapping
STEP=$((STEP+1))

# run bowtie mapping to mRNA
# inside the script, only >22mer are used
BOWTIE_INDEX=${COMMON_FOLDER}/mRNA/mRNA
[ ! -f .${JOBUID}.status.${STEP}.mRNA ] && \
	mkdir -p mRNA && \
	direct_mapping3.sh \
		${BOWTIE_INDEX} \
		$READS_DIR/${INSERT} \
		${mRNA_MM} \
		mRNA \
		"mRNA" \
		$CPU \
		$NormScale && \
	touch .${JOBUID}.status.${STEP}.mRNA
STEP=$((STEP+1))

################
# Joining Pdfs #
################
## by default, only reporting the piRNA feature
echo -e "`date "+$ISO_8601"`\tconcatenating pdfs" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.merge_pdfs ] && \
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX%}.all.pdf \
		$PDF_DIR/`basename ${GENOME_UNIQUEMAP_BED2%.bed*}`.+hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_UNIQUEMAP_BED2%.bed*}`.x_hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_ALLMAP_BED2%.bed*}`.+hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_ALLMAP_BED2%.bed*}`.x_hairpin.lendis.pdf \
		${GENOME_ALLMAP_BED2%bed2}L23mer.features.pdf && \
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX%}.piRNA.pdf \
		$PDF_DIR/`basename ${GENOME_UNIQUEMAP_BED2%.bed*}`.+hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_UNIQUEMAP_BED2%.bed*}`.x_hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_ALLMAP_BED2%.bed*}`.+hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_ALLMAP_BED2%.bed*}`.x_hairpin.lendis.pdf \
		${GENOME_ALLMAP_BED2%bed2}L23mer.features.pdf && \
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX%}.siRNA.pdf \
		$PDF_DIR/`basename ${GENOME_UNIQUEMAP_BED2%.bed*}`.+hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_UNIQUEMAP_BED2%.bed*}`.x_hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_ALLMAP_BED2%.bed*}`.+hairpin.lendis.pdf \
		$PDF_DIR/`basename ${GENOME_ALLMAP_BED2%.bed*}`.x_hairpin.lendis.pdf \
		${GENOME_ALLMAP_BED2%bed2}20-22mer.features.pdf && \
	touch  .${JOBUID}.status.${STEP}.merge_pdfs
STEP=$((STEP+1))

###############
# Cleaning up #
###############
#gzip -f $GENOME_ALLMAP_BED2
#gzip -f $GSNAP_ALLMAP_BED12
#gzip -f $GENOME_UNIQUEMAP_BED2
#gzip -f $GENOME_UNIQUEMAP_BED2_hairpin
#gzip -f ${GENOME_ALLMAP_BED2%bed2}L23mer.bed2
#gzip -f ${GENOME_ALLMAP_BED2%bed2}20-22mer.bed2
#gzip -f ${GENOME_ALLMAP_BED2%bed2}l20mer.bed2

# count in-cluster/out-cluster trn contribution to mature piRNA
case ${ORGANISM} in
fly)
	mkdir -p OSC_Masked
# temp solution
	FLY_SOMA_MASKED_INDEX=/home/hanb/scratch/cd/smallRNA_pipeline_output2/pipeline_output/Brennecke.SRA.WT.unox.OSC/mask/dm3.OSC_masked
	# remapping
	awk 'length($1)>22' ${INPUT%.insert}.${Genome}v${genome_MM}a.al.insert > ${INPUT%.insert}.${Genome}v${genome_MM}a.al.L22.insert && \
	bowtie -S -p $CPU -r --best --strata -a -v 0 $FLY_SOMA_MASKED_INDEX ${INPUT%.insert}.${Genome}v${genome_MM}a.al.L22.insert | samtools view -bSF 0x4 - 2>/dev/null | bedtools bamtobed -i - > ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.bed && \
	istBedToBed2 ${INPUT%.insert}.${Genome}v${genome_MM}a.al.L22.insert ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.bed > ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.bed2 && \
	rm -rf ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.bed
	#FLY_PIRNA_CLUSTER=/home/hanb/nearline/small_RNA_Pipeline/common_files/dm3/UCSC_BEDS/Brennecke.pirnaCluster.bed
	FLY_PIRNA_CLUSTER=/home/wangw1/pipeline_dm3/common/Brennecke.pirnaCluster.x2.xflam.bed # exclude flam and cluster2 to reduce the effect of soma
	#FLY_TRN_IN_CLUSTER=/home/hanb/nearline/small_RNA_Pipeline/common_files/dm3/UCSC_BEDS/transposon.inCluster.bed2
	FLY_TRN_IN_CLUSTER=/home/wangw1/pipeline_dm3/common/transposon.group.bed6.incluster.xflam.xcluster2
	FLY_TRN_OUT_CLUSTER=/home/hanb/nearline/small_RNA_Pipeline/common_files/dm3/UCSC_BEDS/transposon.outCluster.bed2
	siRNA_SUM_TABLE=$SUMMARY_DIR/`basename ${GENOME_ALLMAP_BED2%.bed2}.siRNA.summary`
	siRNA_DEPTH=$(awk '{if ($1=="FLY_cisNATs_all_reads"||$1=="FLY_STRUCTURE_LOCI_all_reads") a+=$2}END{print a/1000000}' $siRNA_SUM_TABLE)
	
	echo -e "`date "+$ISO_8601"`\tin-cluster; out-cluster analysis" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.in_cluster_trn_count ] && \
	bedtoolswo intersect -wa -u -a ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.bed2 -b $FLY_PIRNA_CLUSTER > ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_cluster.bedwo && \
	awk 'BEGIN{FS=OFS="\t"}{if (ARGIND==1) { ct[$7]++ } else { if (ct[$7]==$5) print $1,$2,$3,$4,$5,$6; } }' ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_cluster.bedwo  ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_cluster.bedwo | tee ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.bed | \
	bedtoolswo intersect -wo -a stdin -b $FLY_TRN_IN_CLUSTER | tee ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.bedwo | \
	awk -v depth=$siRNA_DEPTH '{if ($6==$(NF-2)) a[$10]+=$(NF-1)*$4/$NF/$5; else b[$10]+=$(NF-1)*$4/$NF/$5; total[$10]=1}END{for (c in total) { printf "%s\t%.2f\t%.2f\n", c, (a[c]?a[c]:0)/depth, (b[c]?b[c]:0)/depth}}' | \
	sort -k1,1 | \
	tee ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.trn_count | \
	awk '{split ($1,a,"."); p[a[1]]+=$2; n[a[1]]+=$3}END{for (t in p) printf "%s\t%.2f\t%.2f\n", t, p[t], n[t]}' | \
	sort -k1,1 \
	> ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.trn_count.pooled && \
	Rscript $PIPELINE_DIRECTORY/bin/scatter_plot_trn_contribution.R  ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.trn_count.pooled  ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.trn_contribution && \
	awk -v depth=$siRNA_DEPTH '{if ($6==$(NF-2)) a[$16]+=$(NF-1)*$4/$NF/$5; else b[$16]+=$(NF-1)*$4/$NF/$5; total[$16]=1}END{for (c in total) { printf "%s\t%.2f\t%.2f\n", c, (a[c]?a[c]:0)/depth, (b[c]?b[c]:0)/depth}}' ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.bedwo \
	> ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.cluster_count.pooled && \
	Rscript $PIPELINE_DIRECTORY/bin/scatter_plot_cluster_contribution.R  ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.cluster_count.pooled ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_in_cluster_unique.cluster_contribution && \
	touch .${JOBUID}.status.${STEP}.in_cluster_trn_count $PREFIX 
	
[ ! -f .${JOBUID}.status.${STEP}.out_cluster_trn_count ] && \
	bedtoolswo intersect -v -a ${GENOME_ALLMAP_BED2%bed2}L23mer.bed2 -b $FLY_PIRNA_CLUSTER > ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_cluster.bedv && \
	awk 'BEGIN{FS=OFS="\t"}{if (ARGIND==1) { ct[$7]++ } else { if (ct[$7]==$5) print $1,$2,$3,$4,$5,$6; } }' ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_cluster.bedv  ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_cluster.bedv | tee ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_out_cluster_unique.bed | \
	bedtoolswo intersect -wo -a stdin -b $FLY_TRN_OUT_CLUSTER | \
	awk -v depth=$siRNA_DEPTH '{if ($6==$(NF-2)) a[$10]+=$(NF-1)*$4/$NF/$5; else b[$10]+=$(NF-1)*$4/$NF/$5; total[$10]=1}END{for (c in total) { printf "%s\t%.2f\t%.2f\n", c, (a[c]?a[c]:0)/depth, (b[c]?b[c]:0)/depth}}' | \
	sort -k1,1 | \
	tee ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_out_cluster_unique.trn_count | \
	awk '{split ($1,a,"."); p[a[1]]+=$2; n[a[1]]+=$3}END{for (t in p) printf "%s\t%.2f\t%.2f\n", t, p[t], n[t]}' | \
	sort -k1,1 \
	> ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_out_cluster_unique.trn_count.pooled && \
	Rscript $PIPELINE_DIRECTORY/bin/scatter_plot_trn_contribution.R  ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_out_cluster_unique.trn_count.pooled  ${OUTDIR}/OSC_Masked/${PREFIX}.masked.v0a.piRNA_out_cluster_unique.trn_contribution && \
	touch .${JOBUID}.status.${STEP}.out_cluster_trn_count
;;
esac

# print out marker for finishing pipeline
rm -rf ${OUTDIR}"/.finished_*" 
touch  ${OUTDIR}"/.finished_${Genome}_${smallRNA_Pipeline_Version}"
echo -e "`date "+$ISO_8601"`\trun finished" | tee -a $LOG
echo -e "---------------------------------" | tee -a $LOG


