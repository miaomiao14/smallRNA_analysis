#!/bin/bash -x

##smRNA Pipeline for Bombyx mori version 1
##Download the annotation data from http://www.silkdb.org/silkdb/doc/download.html
##10-23-2013
##Author: Wei Wang(wei.wang2@umassmed.edu)


#silkworm genome v2 bowtie index
#/home/wangw1/src/bowtie-0.12.9/indexes/silkwormv2



#1. map to the rRNAs, tRNAs and other nncRNAs

#2.map to the genome 
echo "mapping to the genome...."  >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
/home/wangw1/git/smallRNA_analysis/bm_smRNA_pipeline/run_bowtie_bm.pl $inserts.uniq.reads $mm /home/wangw1/src/bowtie-0.12.9/indexes/silkwormv2 match2_all.out
uniqmap.pl $inserts.match2_all.out > $inserts.uniqmap.match2_all.out
echo `date` >> $dir/output/$insertsname.log
echo "mapping to the genome done" >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log

cut -f1,2 $inserts.match2_all.out  |uniq.lines+ 0 > $inserts.match2_all.out.uniq.reads
cut -f1,2 $inserts.uniqmap.match2_all.out  |uniq.lines+ 0 > $inserts.uniqmap.match2_all.out.uniq.reads

#3. intersect with transposons,genes,miRNAs,









################
# Major Config #
################
# pipeline version
export smallRNA_Pipeline_Version=1.0.0
# this pipeline is still in debug mode
export DEBUG=1 
# pipeline address: if you copy all the files to another directory, this is the place to change; under this directory sits two directories, bin and common. bin stores all the binary executables and common stores all the information of each ORGANISM.
export PIPELINE_DIRECTORY=/home/hanb/nearline/small_RNA_Pipeline
# set PATH to be aware of pipeline/bin; this bin directory needs to be searched first
export PATH=${PIPELINE_DIRECTORY}/bin:$PATH

#########
# USAGE #
#########
# usage function
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -i input_file.fq[.gz] -s fly|mouse -o output_directory[current directory] -c cpu[8] 

This is a single library small RNA pipeline developed in the Zamore Lab in 
University of Massachusetts Medical School. 
Please email Bo.Han@umassmed.edu for any questions or bugs. 
Thanks for using it. 

OPTIONS:
	-h      Show this message
	-i      Input file in fastq or gzipped fastq format, with full directory
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
			FQ=$OPTARG
		;;
		o)
			OUTDIR=$OPTARG
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
checkExist "genomeCoverageBed2" # modified bedtools to take in-house bed2 format and make bedgraph format
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
echo -e "\e[1;35mDone with testing required softwares/scripts, starting pipeline...\e[0m"

##################
# File/Dir check #
##################
# check whether input file exist
[ ! -f $FQ ] && echo -e "\e[1;31mError: Cannot find input file $FQ\e[0m" && exit 1
# if CPU is undefined or containing non-numeric char, then use 8
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
# get full path of the FQ and use the full path when getting the insert file.
FULLPATH_FQ=`readlink -f $FQ`
# use $FQ later only as filename
FQ=${FQ##*/}
# OUTPUT directory defined? otherwise use current directory
[ ! -z $OUTDIR ] || OUTDIR=$PWD
# test wether we need to create the new directory
mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the direcory of input fastq file\e[0m"
# enter destination direcotry, since later it will be ran in this directory, there is NO need to add OUTPUT DIR in the following codes
cd ${OUTDIR} || (echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" && exit 1)
# test writtability
touch .writting_permission && rm -rf .writting_permission || (echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" && exit 1)

#############
# Variables #
#############
# step counter
STEP=1
# unique job id
JOBUID=`echo ${FQ##*/} | md5sum | cut -d" " -f1`
# formating date for log
ISO_8601='%Y-%m-%d %H:%M:%S %Z'
# insert file, a format with two fields delimited by a tab. Sequence and number of times it was read, used to save time/space; quality information is lost
INSERT=${FQ%.f[qa]*}.insert
# table to store the basic statistics of the library (genomic mappability). PS: counts of each genomic features will be stored in .summary file, produced by intersect_all.sh 
TABLE=${FQ%.f[qa]*}.basic_stats
# log file to store the timing
LOG=${FQ%.f[qa]*}.pipeline_${smallRNA_Pipeline_Version}.log
# set ORGANISM, convert upper cases into lower case so it is easier to write the case statement later
ORGANISM=`echo ${ORGANISM} | tr '[A-Z]' '[a-z]'`
# assign different values to the generalized variables (same name for different ORGANISMs) according to which ORGANISM fed
case ${ORGANISM} in
human) #TODO
	Genome=hg19
	COMMON_FOLDER=$PIPELINE_DIRECTORY/common_files/$Genome
	declare -a piRNACLUSTER=(" ")
	declare -a transposon_list=(" ")
	rRNA_MM=2
	hairpin_MM=1
	genome_MM=1
	cluster_MM=1
	transposon_MM=3
;;
mouse)
# using mm9 for mouse 
	Genome=mm9
	COMMON_FOLDER=$PIPELINE_DIRECTORY/common_files/$Genome
## mouse piRNA cluster (bowtie index)
	piRNA_cluster_pachytene=$COMMON_FOLDER/piRNA.cluster.pachytene
	piRNA_cluster_hybrid=$COMMON_FOLDER/piRNA.cluster.hybrid
	piRNA_cluster_prepachytene=$COMMON_FOLDER/piRNA.cluster.prepachytene
# we use an array here so that when we add or delete any members, we only have to change one time
	declare -a piRNACLUSTER=("piRNA_cluster_prepachytene" "piRNA_cluster_hybrid" "piRNA_cluster_pachytene")
# mouse transposons to be checked individually (bowtie index)
	L1=$COMMON_FOLDER/L1
	IAP=$COMMON_FOLDER/IAP
	declare -a transposon_list=("L1" "IAP")
# mismatch configuration variables (for bowtie and gsnap)
	rRNA_MM=3
	hairpin_MM=1
	genome_MM=1
	cluster_MM=1
	transposon_MM=3
;;
mm10)
	Genome=mm10
	COMMON_FOLDER=$PIPELINE_DIRECTORY/common_files/$Genome
	declare -a piRNACLUSTER=(" ")
	declare -a transposon_list=(" ")
	rRNA_MM=3
	hairpin_MM=1
	genome_MM=1
	cluster_MM=1
	transposon_MM=3
;;
fly)
# using dm3 for fly
	Genome=dm3
	COMMON_FOLDER=$PIPELINE_DIRECTORY/common_files/$Genome
# fly piRNA cluster (bowtie index)
	piRNA_cluster=$COMMON_FOLDER/piRNA_cluster
	declare -a piRNACLUSTER=("piRNA_cluster")
# fly transposon (bowtie index)
	gypsy=$COMMON_FOLDER/gypsy
	gypsy12=$COMMON_FOLDER/gypsy12
	HetA=$COMMON_FOLDER/Het-A
	declare -a transposon_list=("gypsy" "gypsy12" "HetA")
# mismatch configuration variables
	rRNA_MM=3
	hairpin_MM=0
	genome_MM=0
	cluster_MM=1
	transposon_MM=2
;;
bombyx)
# using dm3 for fly
	Genome=bmv2
	COMMON_FOLDER=/home/wangw1/pipeline_bm/common
# Bombyx piRNA cluster (bowtie index)
	#piRNA_cluster=$COMMON_FOLDER/piRNA_cluster
	#declare -a piRNACLUSTER=("piRNA_cluster")
# Bombyx transposon (bowtie index)
	#gypsy=$COMMON_FOLDER/gypsy
	#gypsy12=$COMMON_FOLDER/gypsy12
	#HetA=$COMMON_FOLDER/Het-A
	#declare -a transposon_list=("gypsy" "gypsy12" "HetA")
# mismatch configuration variables
	rRNA_MM=3
	hairpin_MM=0
	genome_MM=0
	cluster_MM=1
	transposon_MM=2
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
# premappinglist mapping happens after rRNA and miRNA hairpin mapping (is it appropriate?)
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
mm10)
	declare -a PreMappingList=(" ") 
	declare -a PreMappingMM=(" ")
;;
fly)
	dmel_virus=$COMMON_FOLDER/dmel_virus
	declare -a PreMappingList=("dmel_virus") 
	declare -a PreMappingMM=("2")
;;
bombyx)

	declare -a PreMappingList=("") 
	declare -a PreMappingMM=("")	
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
	awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}'  $1 | sort -k1,1n 
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
[ ! -f .${JOBUID}.status.${STEP}.rRNA_mapping ] && \
	bowtie -r -S -v $rRNA_MM -k 1 -p $CPU \
		--un ${FQ%.f[qa]*}.x_rRNA.insert \
		$BOWTIEINDEX_rRNA \
		${INSERT} \
		1> /dev/stdout \
		2> /dev/null | \
	samtools view -bSF 0x4 - 2>/dev/null | \
	bedtools bamtobed -i - > ${FQ%.f[qa]*}.rRNA.insertBed && \
	istBedToBed2 ${INSERT} ${FQ%.f[qa]*}.rRNA.insertBed > ${FQ%.f[qa]*}.rRNA.bed2 && rm -rf ${FQ%.f[qa]*}.rRNA.insertBed && \
	totalReads=`awk '{a+=$2}END{printf "%d", a}' ${INSERT}` && \
	rRNAReads=`bedwc ${FQ%.f[qa]*}.rRNA.bed2` && \
	nonrRNAReads=$((totalReads-rRNAReads)) && \
    touch .${JOBUID}.status.${STEP}.rRNA_mapping # by checking all the status, ll -a | grep status, you will know which step has finished 
## increment the step counter, no matter previous step succeed or not. Using a variable will make it easy to insert/delete steps
STEP=$((STEP+1))
    
#########################
# miRNA hairpin Mapping #
#########################
# mapping to miRNA hairpins
echo -e "`date "+$ISO_8601"`\tmapping to hairpin, with $hairpin_MM mismatch(es) allowed" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.hairpin_mapping ] && \
	bowtie -r -v $hairpin_MM -a --best --strata -p $CPU -S \
		${BOWTIE_PHRED_OPTION} \
		--al ${FQ%.f[qa]*}.x_rRNA.hairpin.insert \
		--un ${FQ%.f[qa]*}.x_rRNA.x_Hairpin.insert \
		$BOWTIEINDEX_Hairpin \
		${FQ%.f[qa]*}.x_rRNA.insert \
		1> /dev/stdout \
		2> /dev/null  | \
	samtools view -bSF 0x4 - 2>/dev/null | \
	bedtools bamtobed -i - | awk '$6=="+"' > ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.bed && \
	istBedToBed2 ${FQ%.f[qa]*}.x_rRNA.insert ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.bed > ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.bed2 && \
	rm -rf ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.bed && \
	hairpinReads=`bedwc ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.bed2` && \
	bed2lendis ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.bed2 > ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.lendis && \
	bowtie -r -v $genome_MM -a --best --strata -p $CPU \
		-S \
		$BOWTIEINDEX_Genome \
		${FQ%.f[qa]*}.x_rRNA.hairpin.insert \
		1> /dev/stdout \
		2> ${FQ%.f[qa]*}.x_rRNA.hairpin.${Genome}v${genome_MM}m1.log | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	bedtools bamtobed -i - > ${FQ%.f[qa]*}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed && \
	istBedToBed2 ${FQ%.f[qa]*}.x_rRNA.hairpin.insert ${FQ%.f[qa]*}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed > ${FQ%.f[qa]*}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed2 && \
	rm -rf ${FQ%.f[qa]*}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed && \
	touch .${JOBUID}.status.${STEP}.hairpin_mapping
STEP=$((STEP+1))
# run miRNA pipeline
echo -e "`date "+$ISO_8601"`\trunning miRNA analysis pipeline" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.miRNA_pipeline ] && \
	miRNA_pipeline.sh $BOWTIEINDEX_Hairpin $BOWTIEINDEX_Mature ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.bed2  && \
	touch .${JOBUID}.status.${STEP}.miRNA_pipeline
STEP=$((STEP+1))
	
######################
# Pre Genome Mapping #
######################
# count how many indexed need to map
COUNTER=0
# running variable
INPUT=${FQ%.f[qa]*}.x_rRNA.x_Hairpin.insert
# if there are any indexes need to be run
[[ ${#PreMappingList[@]} > 0 ]] && \
for COUNTER in `seq 0 $((${#PreMappingList[@]}-1))`; do \
	TARGET=${PreMappingList[$COUNTER]} && \
	MM=${PreMappingMM[$COUNTER]} && \
	echo -e "`date "+$ISO_8601"`\tmapping to ${TARGET}, with $MM mismatch(es) allowed" | tee -a $LOG && \
	[ ! -f .${JOBUID}.status.${STEP}.${TARGET}_mapping ] && \
		bowtie -r -v $MM -a --best --strata -p $CPU -S \
			${BOWTIE_PHRED_OPTION} \
			--un ${INPUT%.f[qa]*}.x_${TARGET}.insert \
			${!TARGET} \
			$INPUT \
			1> /dev/stdout \
			2> /dev/null | \
		samtools view -bSF 0x4 - 2>/dev/null | bedtools bamtobed -i - > ${INPUT%.f[qa]*}.${TARGET}.v${hairpin_MM}a.bed && \
		istBedToBed2 $INPUT ${INPUT%.f[qa]*}.${TARGET}.v${hairpin_MM}a.bed > ${INPUT%.f[qa]*}.${TARGET}.v${hairpin_MM}a.bed2 && \
		rm -rf ${INPUT%.f[qa]*}.${TARGET}.v${hairpin_MM}a.bed && \
		mappingReads=`bedwc ${INPUT%.f[qa]*}.${TARGET}.v${hairpin_MM}a.bed2` && \
		touch .${JOBUID}.status.${STEP}.${TARGET}_mapping
	STEP=$((STEP+1))
# set OUPUT of this loop to the INPUT of next
	INPUT=${INPUT%.f[qa]*}.x_${TARGET}.insert
done

##################
# Genome Mapping #
##################
# take the OUTPUT of last step as INPUT
INSERT=${INPUT}
# bed2 format storing all mappers for genomic mapping
allBed2=${INPUT%.insert}.${Genome}v${genome_MM}.all.bed2
# bed2 format storing unique mappers for genomic mapping
uniqueBed2=${INPUT%.insert}.${Genome}v${genome_MM}.unique.bed2
# bed2 format storing multiple mappers for genomic mapping
multipBed2=${INPUT%.insert}.${Genome}v${genome_MM}.multip.bed2
# bed2 format storing unique mappers for genomic mapping and miRNA hairpin mapper
uniqueBed2_hairpin=${INPUT%.insert}.${Genome}v${genome_MM}.unique.+hairpin.bed2
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
		2> ${INPUT%.insert}.${Genome}v${genome_MM}m1.insert.log | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	bedtools bamtobed -i - > ${INPUT%.insert}.${Genome}v${genome_MM}a.insert.bed && \
	istBedToBed2 $INPUT ${INPUT%.insert}.${Genome}v${genome_MM}a.insert.bed > ${allBed2} && \
	rm -rf ${INPUT%.insert}.${Genome}v${genome_MM}a.insert.bed && \
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
		--query-unk-mismatch=1 \
		--genome-unk-mismatch=1 \
		--maxsearch=1000000 \
		-n 1000000 \
		--quiet-if-excessive \
		-t $CPU \
		-N 0 \
		-s ${COMMON_FOLDER}/${Genome}/${Genome}.splicesites.iit \
		${INPUT%.insert}.${Genome}v${genome_MM}a.un.fa \
		2> ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.log | \
	samtools view -bS -F0x4 - 2>/dev/null | bedtools bamtobed -bed12 -i - > ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed1 && \
	spliceMapTotal=`awk '{split($4,arr,"_"); ct[$4]=arr[2]}END{for (a in ct){total+=ct[a]} print total}' ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed1` && \
	awk 'BEGIN{OFS="\t"}NR==FNR{++ct[$4]} {$5=ct[$4]; split($4,arr,"_"); $4=arr[2];print}' ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed1 ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed1 > ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed2 && \
	awk '$5==1' ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed2 | bedtools bed12tobed6 -i stdin > ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.unique.bed6 && \
	rm -rf ${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.bed1 ${INPUT%.insert}.${Genome}v${genome_MM}a.un.fa && \
	touch .${JOBUID}.status.${STEP}.gsnap_genome_mapping
STEP=$((STEP+1))
	
# separating unique and multiple mappers
echo -e "`date "+$ISO_8601"`\tseparating unique and multiple mappers" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.separate_unique_and_multiple ] && \
	awk 'BEGIN{OFS="\t"}{if ($5==1) {print $0 >> "/dev/stdout"} else {print $0 >> "/dev/stderr"}}' ${allBed2} \
		1> ${uniqueBed2} \
		2> ${multipBed2} && \
	uniqueMapCount=`bedwc ${uniqueBed2}` && \
	multipMapCount=`bedwc ${multipBed2}` && \
	touch .${JOBUID}.status.${STEP}.separate_unique_and_multiple
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
		echo "awk '\$3-\$2==19' 				${allBed2} > ${allBed2%bed2}19mer.bed2"			>  $paraFile && \
		echo "awk '\$3-\$2 >19 && \$3-\$2<24'	${allBed2} > ${allBed2%bed2}20-23mer.bed2" 		>> $paraFile && \
		echo "awk '\$3-\$2 >23'					${allBed2} > ${allBed2%bed2}L23mer.bed2"		>> $paraFile && \
		ParaFly -c $paraFile -CPU $CPU && \
		rm -rf ${paraFile}* && \
		touch  .${JOBUID}.status.${STEP}.sep_length
	  STEP=$((STEP+1))
;;
## fly specific analysis
fly)
	echo -e "`date "+$ISO_8601"`\tseparating reads based on length" | tee -a $LOG
	[ ! -f .${JOBUID}.status.${STEP}.sep_length ] && \
		paraFile=${RANDOM}${RANDOM}.para && \
		echo "awk '\$3-\$2==21' ${allBed2} > ${allBed2%bed2}21mer.bed2"   > $paraFile && \
		echo "awk '\$3-\$2 >22' ${allBed2} > ${allBed2%bed2}L22mer.bed2" >> $paraFile && \
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
	awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}' ${allBed2}  | sort -k1,1n > ${allBed2}.lendis && \
	awk '{a[$7]=$4}END{m=0; for (b in a){c[length(b)]+=a[b]; if (length(b)>m) m=length(b)} for (d=1;d<=m;++d) {print d"\t"(c[d]?c[d]:0)}}' ${uniqueBed2}  | sort -k1,1n > ${uniqueBed2}.lendis && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_lendis.R ${allBed2}.lendis ${allBed2%.bed*}.x_hairpin && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_lendis.R ${uniqueBed2}.lendis ${uniqueBed2%.bed*}.x_hairpin && \
	awk '{ct[$1]+=$2}END{for (l in ct) {print l"\t"ct[l]}}' ${allBed2}.lendis    ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.lendis | sort -k1,1n > ${allBed2}.+hairpin.lendis && \
	awk '{ct[$1]+=$2}END{for (l in ct) {print l"\t"ct[l]}}' ${uniqueBed2}.lendis ${FQ%.f[qa]*}.x_rRNA.hairpin.v${hairpin_MM}a.lendis | sort -k1,1n > ${uniqueBed2}.+hairpin.lendis && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_lendis.R ${allBed2}.+hairpin.lendis    ${allBed2%.bed*}.+hairpin && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_lendis.R ${uniqueBed2}.+hairpin.lendis ${uniqueBed2%.bed*}.+hairpin && \
	touch .${JOBUID}.status.${STEP}.plotting_length_dis
STEP=$((STEP+1))


##################
# Print to table #
##################

echo -e "total reads as input of the pipeline\t${totalReads}" > $TABLE && \
echo -e "rRNA reads with ${rRNA_MM} mismatches\t${rRNAReads}" >> $TABLE && \
echo -e "genome mapping reads (-rRNA; +miRNA_hairpin)\t$((uniqueMapCount+multipMapCount+hairpinReads+spliceMapTotal))" >> $TABLE && \
echo -e "miRNA hairpin reads\t${hairpinReads}" >> $TABLE && \
echo -e "genome mapping reads (-rRNA; -miRNA_hairpin)\t$((uniqueMapCount+multipMapCount+spliceMapTotal))" >> $TABLE && \
echo -e "genome mapping reads (-rRNA; -miRNA_hairpin; -exon_exon_junction_mapper)\t$((uniqueMapCount+multipMapCount))" >> $TABLE && \
echo -e "genome unique mapping reads (-rRNA; -miRNA_hairpi; -exon_exon_junction_mapper)\t${uniqueMapCount}" >> $TABLE && \
echo -e "genome multiple mapping reads (-rRNA; -miRNA_hairpi; -exon_exon_junction_mapper)\t${multipMapCount}" >> $TABLE && \
echo -e "exon-exon junction reads (unmappable to genome)\t${spliceMapTotal}" >> $TABLE && \

###############################
# Intersecting Genome Feature #
###############################
# running intersecting bed
echo -e "`date "+$ISO_8601"`\trunning intersectBed for all features" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.intersect_all ] && \
intersect_all.sh \
	${uniqueBed2} \
	${multipBed2} \
	${allBed2%.bed2}.summary \
	$ORGANISM \
	$GENOME_FA \
	$CPU && \
	touch .${JOBUID}.status.${STEP}.intersect_all
STEP=$((STEP+1))

########################
# Makeing BigWig Files #
########################
# normalization factor, currently using unique genome mappers and all hairpin mappers
NormScale=`echo $((uniqueMapCount+hairpinReads+spliceMapTotal)) | awk '{printf "%f",1000000.0/$1}'`
# make BW files
echo -e "`date "+$ISO_8601"`\tmaking bigWig file for genome browser" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.make_bigWig ] && \
cat ${FQ%.f[qa]*}.x_rRNA.hairpin.${Genome}v${genome_MM}a.bed2 ${uniqueBed2} > $uniqueBed2_hairpin && \
	bed22bw.sh \
		${uniqueBed2_hairpin} \
		$CHROM \
		${NormScale} \
		${INPUT%.insert}.${Genome}v${genome_MM}a.un.gsnap.v${genome_MM}.unique.bed6 \
		$CPU && \
	rm -rf $uniqueBed2_hairpin && \
	touch .${JOBUID}.status.${STEP}.make_bigWig
STEP=$((STEP+1))

###################################
# Mapping to transposons directly #
###################################
# run bowtie mapping to some transposons directly
for transposon in ${transposon_list[@]}; do 
	echo -e "`date "+$ISO_8601"`\trunning direct mapping to $transposon, with ${transposon_MM} mismatch(es) allowed " | tee -a $LOG
	[ ! -f .${JOBUID}.status.${STEP}.${transposon} ] && \
		bowtie -r -v ${transposon_MM} -a --best --strata -p $CPU \
			--al ${INSERT%.insert}.${transposon}.v${transposon_MM}a.al.insert -S \
			${!transposon} \
			${INSERT} \
			1> /dev/stdout \
			2> /dev/null | \
		samtools view -uS -F0x4 - 2>/dev/null | bedtools bamtobed -i - > ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.bed  && \
		istBedToBed2 $INSERT ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.bed > ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.bed2 && \
		rm -rf ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.bed && \
		MapReads=`bedwc ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.bed2` && \
		echo -e "${transposon} direct mapper reads (independent of genome mapping)\t${MapReads}" >> $TABLE && \
		ppbed2 -a  ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.bed2 -b ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.bed2 > ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.ppbed2 && \
		Rscript --slave ${PIPELINE_DIRECTORY}/bin/draw_pp.R ${INSERT%.insert}.${transposon}.v${transposon_MM}a.insert.ppbed2 ${INSERT%.insert}.${transposon}.v${transposon_MM}a && \
		touch .${JOBUID}.status.${STEP}.${transposon}
	STEP=$((STEP+1))
 done

#####################################
# Mapping to piRNA cluster directly #
#####################################
# mapping to piRNA cluster directly
for cluster in ${piRNACLUSTER[@]}; do 
	echo -e "`date "+$ISO_8601"`\trunning direct mapping to $cluster, with $cluster_MM mismatch(es) allowed" | tee -a $LOG
	[ ! -f .${JOBUID}.status.${STEP}.${cluster} ] && \
		bowtie -r -v $cluster_MM -a --best --strata -p $CPU \
			-S \
			${!cluster} \
			${INSERT} \
			1> /dev/stdout \
			2> /dev/null | \
		samtools view -uS -F0x4 - 2>/dev/null | bedtools bamtobed -i - > ${INSERT%.insert}.${cluster}.v${cluster_MM}a.insert.bed && \
		istBedToBed2 ${INSERT} ${INSERT%.insert}.${cluster}.v${cluster_MM}a.insert.bed > ${INSERT%.insert}.${cluster}.v${cluster_MM}a.insert.bed2 && \
		rm -rf ${INSERT%.insert}.${cluster}.v${cluster_MM}a.insert.bed && \
		ppbed2each -o /dev/stdout -a ${INSERT%.insert}.${cluster}.v${cluster_MM}a.insert.bed2 -b ${INSERT%.insert}.${cluster}.v${cluster_MM}a.insert.bed2 | \
		awk '{if (ARGIND==1) {a[$4]=1} else {b[$1]=$2}}END{for (c in a) {print c"\t"(b[c]?b[c]:0)}}' ${!cluster}".bed" - > ${INSERT%.insert}.${cluster}.zscores && \
		touch  .${JOBUID}.status.${STEP}.${cluster}
	STEP=$((STEP+1))
 done

################
# Joining Pdfs #
################
echo -e "`date "+$ISO_8601"`\tconcatenating pdfs" | tee -a $LOG
[ ! -f .${JOBUID}.status.${STEP}.merge_pdfs ] && \
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${FQ%.fq}.pdf \
		${uniqueBed2%.bed*}.+hairpin.lendis.pdf \
		${uniqueBed2%.bed*}.x_hairpin.lendis.pdf \
		${allBed2%.bed*}.+hairpin.lendis.pdf \
		${allBed2%.bed*}.x_hairpin.lendis.pdf \
		${uniqueBed2}.features.pdf && \
	touch  .${JOBUID}.status.${STEP}.merge_pdfs
STEP=$((STEP+1))

###############
# Final Table #
###############

###############
# Cleaning up #
###############
[ ! -f .${JOBUID}.status.${STEP}.gzip ] && \
	paraFile="clean.para" && rm -rf $paraFile ${paraFile}.completed && \
	for file in `ls | grep '\.bed2$' && ls | grep insert$`;
		do
			echo "gzip $i" >> $paraFile;
	done
	ParaFly -c $paraFile -CPU $CPU && \
touch .${JOBUID}.status.${STEP}.gzip

# print out marker for finishing pipeline
rm -rf ${OUTDIR}"/.finished_*" 
touch  ${OUTDIR}"/.finished_${Genome}_${smallRNA_Pipeline_Version}"
echo -e "`date "+$ISO_8601"`\trun finished" | tee -a $LOG
echo -e "---------------------------------" | tee -a $LOG



