#!/bin/bash -x

##smRNA Pipeline for Bombyx mori version 1
##Download the annotation data from http://www.silkdb.org/silkdb/doc/download.html
##10-23-2013
##Author: Wei Wang(wei.wang2@umassmed.edu)


#silkworm genome v2 bowtie index
#/home/wangw1/src/bowtie-0.12.9/indexes/silkwormv2



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
#checkExist "genomeCoverageBed2" # modified bedtools to take in-house bed2 format and make bedgraph format
checkExist "weblogo" # make weblogo 
checkExist "information_content.py" # calculate information content from fasta
checkExist "ParaFly" # a handy tool in trinity package to wrap a serials independent commands and run them with multiple threads
checkExist "bedGraphToBigWig" # a tool from kent source to convert bedgraph to bigWig format
checkExist "phred_test" # in-house C++ program to find out wether a fastq format is in phred+33 or phred+64 mode
#checkExist "ppbed2" # in-house C++ program to calculate ping-pong score from two bed2 files
#checkExist "ppbed2each" # in-house C++ program to calculate ping-pong score of each chr/transcript from two bed2 file
#checkExist "istBedToBed2" # in-house C++ program to convert bed file, which is made by bowtie-mapping of insert file to bed2 format
#checkExist "calculate_miRNA_heterogeneity" # in-house C++ program to calculate miRNA heterogeneity from bed2 format
#checkExist "fastq2insert" # in-house C++ program to convert compressed/uncompressed fastq file to insert file
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
	RRNA=${COMMON_FOLDER}/silkworm_rRNA.gff
	TRNA=${COMMON_FOLDER}/silkworm_tRNA.gff
	KNOWNTE=${COMMON_FOLDER}/silkworm_Publicknow_TE.bed
	ReASTE=${COMMON_FOLDER}/silkworm_ReAS_TE.bed
	MIRNA=${COMMON_FOLDER}/silkworm_miRNA.bed
	GENE=${COMMON_FOLDER}/silkworm_glean.bed
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
BOWTIEINDEX_Genome=$COMMON_FOLDER/$Genome/silkwormv2
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


STEP=1
########################################################
#remove adapter and convert the fastq format to inserts
########################################################


file=${FQ##*/}
filename=`basename $file .fastq`
echo -e "`date "+$ISO_8601"`\tremove the adaptors and convert the format to inserts " | tee -a $LOG
touch .status.${STEP}.preprocessing
[ ! -f .status.${STEP}.preprocessing ] && \
grep -A 1 "@" ${FULLPATH_FQ} | grep -v "@" | grep -v "\-\-"  | uniq.reads+   >${filename}.raw && \
Extract_insert_10mer.pl ${filename}.raw $ADAPTER > ${OUTDIR}/${filename}.insertsout && \
inserts2uniqreads.pl ${OUTDIR}/${filename}.insertsout 18 30 > ${OUTDIR}/${filename}.inserts.trimmed && \
rm ${OUTDIR}/${filename}.insertsout && \
	touch .status.${STEP}.preprocessing
STEP=$((STEP+1))


##################
# Genome Mapping #
##################
# take the OUTPUT of last step as INPUT
INSERT=${OUTDIR}/${filename}.inserts.trimmed
# bed2 format storing all mappers for genomic mapping
allBed2=${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}.all.bed2
# bed2 format storing unique mappers for genomic mapping
uniqueBed2=${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}.unique.bed2
# bed2 format storing multiple mappers for genomic mapping
multipBed2=${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}.multip.bed2
# bed2 format storing unique mappers for genomic mapping and miRNA hairpin mapper
uniqueBed2_hairpin=${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}.unique.+hairpin.bed2
# mapping insert file to genome
echo -e "`date "+$ISO_8601"`\tmapping to genome, with ${genome_MM} mismatch(es) allowed " | tee -a $LOG
[ ! -f .status.${STEP}.genome_mapping ] && \
	bowtie -r -v $genome_MM -a --best --strata -p $CPU \
		--al  ${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}a.al.insert \
		--un  ${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}a.un.insert \
		-S \
		$BOWTIEINDEX_Genome \
		${INSERT} \
		1> /dev/stdout \
		2> ${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}m1.insert.log | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	bedtools bamtobed -i - > ${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}a.insert.bed && \
	istBedToBed2 $INSERT ${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}a.insert.bed > ${allBed2} && \
	rm -rf ${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}a.insert.bed && \
	touch .status.${STEP}.genome_mapping
STEP=$((STEP+1))

#stat
insertsReadNum=`sumcol ${INSERT} 2`
genomeMapReadNum=`sumcol ${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}a.al.insert 2`

insertsSpeciesNum=`wc -l ${INSERT} |cut -f1 -d" "`
genomeMapSpeciesNum=`wc -l ${INSERT%.inserts.trimmed}.${Genome}v${genome_MM}a.al.insert|cut -f1 -d" "`


#x rRNAs,tRNAs by intersect with rRNA and tRNA annotations,remove both sense and antisense mappers
echo -e "`date "+$ISO_8601"`\tremove rRNAs and tRNAs, remove miRNA mappers" | tee -a $LOG
[ ! -f .status.${STEP}.xrRNA.xtRNA ] && \
[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.bed2 ] && bedtools intersect -a ${allBed2} -b ${RRNA} -v -f 0.99 |bedtools intersect -a - -b ${TRNA} -v -f 0.99 >${allBed2%*.bed2}.all.xrRNA.xtRNA.bed2 && \
touch ${OUTDIR}/.status.${STEP}.xrRNA.xtRNA
STEP=$((STEP+1))


#stat
TRRNA=${COMMON_FOLDER}/silkworm_rRNA_tRNA.gff
[ ! -f ${TRRNA} ] && cat ${RRNA} ${TRNA} >${TRRNA}
[ ! -f ${allBed2%*.bed2}.all.rRNA.tRNA.bed2 ] && bedtools intersect -a ${allBed2} -b ${TRRNA} -wb -f 0.99 >${allBed2%*.bed2}.all.rRNA.tRNA.bed2
[ ! -f ${allBed2%*.bed2}.all.rRNA.tRNA.uniq.reads ] && awk '{OFS="\t"}{print $7,$4}' ${allBed2%*.bed2}.all.rRNA.tRNA.bed2 |sort -u >${allBed2%*.bed2}.all.rRNA.tRNA.uniq.reads

nncReadNum=`sumcol ${allBed2%*.bed2}.all.rRNA.tRNA.uniq.reads 2`
nncSpeciesNum=`wc -l ${allBed2%*.bed2}.all.rRNA.tRNA.uniq.reads|cut -f1 -d" "`

#x miRNAs
echo -e "`date "+$ISO_8601"`\tremove miRNA mappers" | tee -a $LOG
[ ! -f .status.${STEP}.xmiRNA ] && \
[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.bed2 ] && bedtools intersect -a ${allBed2%*.bed2}.all.xrRNA.xtRNA.bed2 -b ${MIRNA} -v -f 0.99 >${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.bed2 
#intersect with miRNAs
[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.S.bed2 ] && bedtools intersect -a ${allBed2%*.bed2}.all.xrRNA.xtRNA.bed2 -b ${MIRNA} -f 0.99 -s -wb >${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.S.bed2 
[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.AS.bed2 ] && bedtools intersect -a ${allBed2%*.bed2}.all.xrRNA.xtRNA.bed2 -b ${MIRNA} -f 0.99 -S -wb >${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.AS.bed2 && \
touch ${OUTDIR}/.status.${STEP}.xmiRNA
STEP=$((STEP+1))


#stat
[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.S.uniq.reads ] && awk '{OFS="\t"}{print $7,$4}'  ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.S.bed2 |sort -u > ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.S.uniq.reads
[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.AS.uniq.reads ] && awk '{OFS="\t"}{print $7,$4}'  ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.AS.bed2 |sort -u > ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.AS.uniq.reads

miRNASReadNum=`sumcol ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.S.uniq.reads 2`
miRNAASReadNum=`sumcol ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.AS.uniq.reads 2`

miRNASSpeciesNum=`wc -l ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.S.uniq.reads|cut -f1 -d" "`
miRNAASSpeciesNum=`wc -l ${allBed2%*.bed2}.all.xrRNA.xtRNA.miRNA.AS.uniq.reads|cut -f1 -d" "`
#length distribution


#interesect with gene, TE annotations
touch ${OUTDIR}/.status.${STEP}.intersectGeneandTE
declare -a TARGETS=("GENE" "KNOWNTE" "ReASTE")
echo -e "`date "+$ISO_8601"`\tintersect with gene and TES, converting format to mapper2" | tee -a $LOG
[ ! -f .status.${STEP}.intersectGeneandTE ] && \
parafly_file=${OUTDIR}/intersectTogeneTE.para && \
for t in ${TARGETS[@]}
do
echo -e "[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.bed2 ] && bedtools intersect -a ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.bed2 -b ${!t} -f 0.99 -s -wb >${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.bed2 " >> $parafly_file ; 
echo -e "[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.bed2 ] && bedtools intersect -a ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.bed2 -b ${!t} -f 0.99 -S -wb >${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.bed2 " >> $parafly_file ;
#echo -ne "cat ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.mapper2 ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.mapper2> ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.mapper2 &&  ">> $parafly_file ;
#echo -e "gzip  ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.mapper2 && rm ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.mapper2 && rm ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.mapper2">> $parafly_file ;
done
if [[ ! -s ${parafly_file}.completed ]] || [[ -f ${parafly_file}.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU 8 -failed_cmds ${parafly_file}.failed_commands
fi
[ $? == 0 ] && \
touch ${OUTDIR}/.status.${STEP}.intersectGeneandTE
STEP=$((STEP+1))

touch ${OUTDIR}/.status.${STEP}.bed2mapper2
[ ! -f .status.${STEP}.bed2mapper2 ] && \
parafly_file=${OUTDIR}/bed2mapper2.para && \
for t in ${TARGETS[@]}
do
echo -e "[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.mapper2 ] && bedtools sort -i ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.bed2  |bedtools groupby -i - -g 1,2,3,4,5,6,7 -c 8,9,10,11,12,13 -o collapse,collapse,collapse,collapse,collapse,collapse | \
			awk 'BEGIN{OFS=\"\\\t\"}{k=split(\$11,kr,\",\");l=split(\$11,lr,\",\"); \$2+=1; for(i=1;i<=k;i++){print \$7,\$4,\$1\":\"\$2\"-\"\$3\"(\"\$6\")\",\"sense\",kr[i],lr[i],\$5,k} }'  \
			> ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.mapper2 ">> $parafly_file ;
echo -e "[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.mapper2 ] && bedtools sort -i ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.bed2  |bedtools groupby -i - -g 1,2,3,4,5,6,7 -c 8,9,10,11,12,13 -o collapse,collapse,collapse,collapse,collapse,collapse | \
			awk 'BEGIN{OFS=\"\\\t\"}{k=split(\$11,kr,\",\");l=split(\$12,lr,\",\"); \$2+=1; for(i=1;i<=k;i++){print \$7,\$4,\$1\":\"\$2\"-\"\$3\"(\"\$6\")\",\"antisense\",kr[i],lr[i],\$5,k} }'  \
			> ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.mapper2 ">> $parafly_file ;
#echo -ne "cat ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.mapper2 ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.mapper2> ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.mapper2 &&  ">> $parafly_file ;
#echo -e "gzip  ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.mapper2 && rm ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.mapper2 && rm ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.mapper2">> $parafly_file ;
done
if [[ ! -s ${parafly_file}.completed ]] || [[ -f ${parafly_file}.failed_commands ]]
then
	ParaFly -c $parafly_file -CPU 8 -failed_cmds ${parafly_file}.failed_commands
fi
[ $? == 0 ] && \
touch ${OUTDIR}/.status.${STEP}.bed2mapper2
STEP=$((STEP+1))

#stat
touch .status.${STEP}.uniqreads.stat
[ ! -f .status.${STEP}.uniqreads.stat ] && \
for t in ${TARGETS[@]}
do
[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.uniq.reads ] && awk '{OFS="\t"}{print $1,$2}' ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.mapper2 |sort -u >${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.uniq.reads
[ ! -f ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.uniq.reads ] && awk '{OFS="\t"}{print $1,$2}' ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.mapper2 |sort -u >${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.uniq.reads
eval SReadNum${t}=`sumcol ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.uniq.reads 2`
eval ASReadNum${t}=`sumcol ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.uniq.reads 2`

eval SSpeciesNum${t}=`wc -l ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.S.uniq.reads|cut -f1 -d" "`
eval ASSpeciesNum${t}=`wc -l ${allBed2%*.bed2}.all.xrRNA.xtRNA.xh.${t}.AS.uniq.reads |cut -f1 -d" "`
done
[ $? == 0 ] && \
#wrote the 
READSTAT=${OUTDIR}/${filename}.reads.stat && \
echo -e "inserts\tgenome_mapping\trRNA&tRNAs\tmiRNA_S\tmiRNA_AS\tGENE_S\tGENE_AS\tKnownTE_S\tKnownTE_AS\tReASTE_S\tReASTE_AS\n" >${READSTAT} && \
echo -e "${insertsReadNum}\t${genomeMapReadNum}\t${nncReadNum}\t${miRNASReadNum}\t${miRNAASReadNum}\t${SReadNumGENE}\t${ASReadNumGENE}\t${SReadNumKNOWNTE}\t${ASReadNumKNOWNTE}\t${SReadNumReASTE}\t${ASReadNumReASTE}\n" >>${READSTAT}	&& \
SPECIESSTAT=${OUTDIR}/${filename}.species.stat && \
echo -e "inserts\tgenome_mapping\trRNA&tRNAs\tmiRNA_S\tmiRNA_AS\tGENE_S\tGENE_AS\tKnownTE_S\tKnownTE_AS\tReASTE_S\tReASTE_AS\n" >${SPECIESSTAT} && \
echo -e "${insertsSpeciesNum}\t${genomeMapSpeciesNum}\t${nncSpeciesNum}\t${miRNASSpeciesNum}\t${miRNAASSpeciesNum}\t${SSpeciesNumGENE}\t${ASSpeciesNumGENE}\t${SSpeciesNumKNOWNTE}\t${ASSpeciesNumKNOWNTE}\t${SSpeciesNumReASTE}\t${ASSpeciesNumReASTE}\n" >>${SPECIESSTAT} && \	

touch ${OUTDIR}/.status.${STEP}.uniqreads.stat
STEP=$((STEP+1))

#seqLogo
touch .status.${STEP}.seqLogo
[ ! -f .status.${STEP}.seqLogo ] && \
for i in `ls *.uniq.reads.gz`
do
ex $i
awk '{OFS="\t"}{if(length($1)>=23 && length($1)<=30) print $0}' ${i%.gz} >${i%.gz}.23-29 &&
seqlogo_ww ${i%.gz}.23-29 29 /home/wangw1/isilon_temp/BmN4/seqLogo r && rm /home/wangw1/isilon_temp/BmN4/seqLogo/*.uniq.seq

done

[ $? == 0 ] && \
touch .status.${STEP}.seqLogo
STEP=$((STEP+1))

[ -f ${allBed2} ] && rm ${allBed2}
for i in `ls *.uniq.reads`
do 
	echo -e "gzip $i" >>${OUTDIR}/parafile.gzip
done
ParaFly -c ${OUTDIR}/parafile.gzip -CPU 8 
gzip ${INSERT}
gzip ${INSERT%.inserts.trimmed}.raw
gzip ${INSERT%.inserts.trimmed}.fastq
gzip ${INSERT%.inserts.trimmed}.*.insert
gzip ${allBed2%*.bed2}.all.xrRNA.xtRNA.bed2
