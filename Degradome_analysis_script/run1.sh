#! /bin/bash -x
# Single-end RNASeq pipeline (CAGE, SE degradome...) in the Zamore Lab
# 2013-06-12
# Bo W Han (bo.han@umassmed.edu)
# Phillip Zamore Lab
# RNA Therapeutics Institute
# HHMI & University of Massachusetts Medical School

################
# Major Config #
################
# VERSION CONTROL
export SE_PIPELINE_VER=1.0
# IN DEBUG MODE
export DEBUG=1 
# PIPELINE DIRECTORY
export SE_PIPELINE_ADD=/home/hanb/nearline/PE_Pipeline/
# RSEM DIRECTORY
RSEM_DIR=/home/hanb/software/trinityrnaseq_r2013-02-25/trinity-plugins/rsem
# SCRIPTS SEARCHING PATH
export PATH=${SE_PIPELINE_ADD}/bin:$RSEM_DIR:${PATH}

#########
# USAGE #
#########
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -i input.fa -s fly|mouse -c 8 -o output_directory 
This is a single end RNA Sequencing pipeline developped in the Zamore Lab in 
University of Massachusetts Medical School. 
Please email Bo.Han@umassmed.edu for any questions or bugs. 
Thanks for using it. 

OPTIONS:
	-h      Show this message
	-i      Path to the intput fastq
	-s      Species [ fly | mouse ]
	-o      Output Directory [ default = current directory ]
	-c      Number of CPUs to use [ default = 8 ]
EOF
echo -en "\e[0m"
}
# initiate strandless
STRAND=0
while getopts "hi:o:s:c:" OPTION
do
	case $OPTION in
		h)
			usage && exit 1
		;;
		i)
			INPUT_FQ=$OPTARG
		;;
		s)
			ORGANISM=$OPTARG
		;;
		o)
			OUTDIR=$OPTARG
		;;
		c)
			CPU=$OPTARG
		;;
		?)
			usage && exit 1
		;;
	esac
done
# if LEFT or RIGHT or ORGANISM is undefined, print out help and exit
if [[ -z $INPUT_FQ ]] || [[ -z $ORGANISM ]]
then
	usage && exit 1
fi
# if CPU is undefined, then use 8
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
# if OUTDIR is undefined, then use PWD
[ ! -z ${OUTDIR} ] || OUTDIR=$PWD
# if unspecified strandless, then use dUTR method

###########################
# running necessity check #
###########################
function checkExist {
	echo -ne "\e[1;32m\"${1}\" is using: \e[0m" && which "$1"
	[[ $? != 0 ]] && echo -e "\e[1;31mError: cannot find software/function ${1}! Please make sure that you have installed the pipeline correctly.\nExiting...\e[0m" && \
	exit 1
}
echo -e "\e[1;35mTesting required softwares/scripts:\e[0m"
checkExist "echo" # used to print message
checkExist "date" # used to get time
checkExist "mkdir" # used to creat new directory
checkExist "rm" # used to delete file/directory
checkExist "mv" # used to move or rename file
checkExist "touch" # used to create new empty file
checkExist "awk" # a programming language
checkExist "samtools" # tools to process sam/bam file
checkExist "bedtools" # tools to process bed formatted file
checkExist "bowtie2" # mapping tool, used for most mapping purpose
checkExist "ParaFly" # a handy tool in trinity package to wrap a serials independent commands and run them with multiple threads
checkExist "intersectPE.sh" # using bedtools bedtopair to assign eash cordant pair to genomic features
checkExist "bed22weblogo.sh" # making information content graph
checkExist "information_content.py" # used by bed22weblogo.sh to get information content
checkExist "weblogo" # draw weblogo
checkExist "bedGraphToBigWig" # convert bedGraph to bigWig
checkExist "bam" # bamtools
echo -e "\e[1;35mDone with testing required softwares/scripts, starting pipeline...\e[0m"

#############
# Variables #
#############
# formating date (in log)
ISO_8601='%Y-%m-%d %H:%M:%S %Z'
# step counter
STEP=1
# set ORGANISM specific files, convert upper cases into lower case so it is easier to write the case statement later
ORGANISM=`echo ${ORGANISM} | tr '[A-Z]' '[a-z]'`
# use the longest PREFIX of $LEFT and $RIGHT as the library name
INPUT_FQ_BASENAME=${INPUT_FQ##*/}
PREFIX=${INPUT_FQ_BASENAME%.f[aq]*}
# Log file
LOG=${PREFIX}.${SE_PIPELINE_VER}.log
# checking read length and determine the best sjdbOverhang
READ_LEN=`head -2 $INPUT_FQ | awk '{getline; printf "%d", length($1)}'`
# overhang
STARsjdbOverhang=$((READ_LEN-1))
# testing phred_format
PHRED_SCORE=`phred_test ${INPUT_FQ}`
case ${PHRED_SCORE} in
33)
	bowtie2PhredOption="--phred33"
	tophatPhredOption=""
;;
64)
	bowtie2PhredOption="--phred64"
	tophatPhredOption="--solexa1.3-quals"
;;
*)
	bowtiePhredOption=""
;;
esac
# decide organism and genomie assembly
case ${ORGANISM} in
mouse)
	Genome=mm9
;;
fly)
	Genome=dm3
;;
*)
	echo -e "\e[1;31mError: Currently only mm9 and dm3 are supported\e[0m"
	usage && exit 1
;;
esac

# other variables
Chrome=${SE_PIPELINE_ADD}/common_files/${Genome}/ChromInfo.txt
# genome fasta
GENOME_FA=${SE_PIPELINE_ADD}/common_files/${Genome}/STARgenomeIndex/genome.fa
# directory with STAR indexes stored
STAR_Genome_INDEX_DIR=${SE_PIPELINE_ADD}/common_files/${Genome}/STARgenomeIndex
# directory with bowtie index for rRNA
STAR_rRNA_INDEX_DIR=${SE_PIPELINE_ADD}/common_files/${Genome}/rRNAIndex
# GTF file // TODO: merge this gtf...
GTF_FILE=${SE_PIPELINE_ADD}/common_files/${Genome}/mouse.testis.00-07.TopHat.Cufflinks.gtf
# BWA index for transposon sequences
BWA_TRANSPOSON_INDEX=${SE_PIPELINE_ADD}/common_files/${Genome}/repBase/transposon.ref
# BWA index for piRNA cluster
BWA_PIRNACLUSTER_INDEX=${SE_PIPELINE_ADD}/common_files/${Genome}/piRNAcluster/piRNAcluster.fa

#######################
#Begin to run pipeline#
#######################
[ ! -f $LOG ] && \
echo -e "`date "+$ISO_8601"`\tbeginning running PE pipeline version $SE_PIPELINE_VER"  | tee -a $LOG || \
echo -e "`date "+$ISO_8601"`\t..resuming running PE pipeline version $SE_PIPELINE_VER" | tee -a $LOG 

# test if STAR index exist, if not, creat it
STAR_Genome_INDEX=$STAR_Genome_INDEX_DIR/$STARsjdbOverhang
[ ! -d $STAR_Genome_INDEX ] && \
	mkdir -p $STAR_Genome_INDEX &&
	STAR --runMode genomeGenerate \
		--runThreadN $CPU \
		--genomeDir $STAR_Genome_INDEX \
		--genomeFastaFiles $GENOME_FA \
		--sjdbOverhang ${STARsjdbOverhang}

########################
#ribosomal rRNA mapping#
########################
echo -e "`date "+$ISO_8601"`\tmapping to rRNA with Bowtie2" | tee -a $LOG
[ ! -f .status.${STEP}.rRNA_mapping ] && \
bowtie2 \
	-x $STAR_rRNA_INDEX_DIR/rRNA \
	-U $INPUT_FQ \
	-q \
	$bowtie2PhredOption \
	--very-fast \
	-k 1 \
	--no-mixed \
	--no-discordant \
	--un ${PREFIX}.x_rRNA.UNmapped.fq \
	-p $CPU \
	-S ${PREFIX}.rRNA.sam \
	2> ${PREFIX}.rRNA.log && \
touch .status.${STEP}.rRNA_mapping
STEP=$((STEP+1))

################
#genome mapping#
################
# set LEFT and RIGHT to be the unmappable reads for rRNA
LEFT=${PREFIX}.x_rRNA.UNmapped.fq
[ ! -f .status.${STEP}.genome_mapping ] && \
STAR \
	--runMode alignReads \
	--genomeDir $STAR_Genome_INDEX \
	--readFilesIn ${LEFT} \
	--runThreadN $CPU \
	--outFilterScoreMin 0 \
	--outFilterScoreMinOverLread 0.72 \
	--outFilterMatchNmin 0 \
	--outFilterMatchNminOverLread 0.72 \
	--outFilterMultimapScoreRange 1 \
	--outFilterMultimapNmax -1 \
	--outFilterMismatchNmax 10 \
	--outFilterMismatchNoverLmax 0.05 \
	--alignIntronMax 0 \
	--alignIntronMin 21 \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--genomeLoad NoSharedMemory \
	--outFileNamePrefix ${PREFIX}.x_rRNA.${Genome}. \
	--outSAMunmapped None \
	--outReadsUnmapped Fastx \
	--outSJfilterReads Unique \
	--seedSearchStartLmax 20 \
	--seedSearchStartLmaxOverLread 1.0 \
	--chimSegmentMin 0 && \
touch .status.${STEP}.genome_mapping
STEP=$((STEP+1))

# getting statistics
InputReads=`grep 'Number of input reads' ${PREFIX}.x_rRNA.${Genome}.Log.final.out | awk '{print $NF}'`
UniquReads=`grep 'Uniquely mapped reads number' ${PREFIX}.x_rRNA.${Genome}.Log.final.out | awk '{print $NF}'`
MultiReads=`grep 'Number of reads mapped to multiple loci' ${PREFIX}.x_rRNA.${Genome}.Log.final.out | awk '{print $NF}'`
AllMapReads=$((UniquReads+MultiReads))
UnMapReads=$((InputReads-UniquReads-MultiReads))
# normalization factor
NormScale=`echo $UniquReads | awk '{printf "%f",1000000.0/$1}'`

# separate unique mapper and multiple mapper
echo -e "`date "+$ISO_8601"`\tprocessing mapped result" | tee -a $LOG
[ ! -f .status.${STEP}.genome_bam_processing ] && \
	samtools view -bS ${PREFIX}.x_rRNA.${Genome}.Aligned.out.sam | tee ${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam | bedtools bamtobed -bed12 -i - > ${PREFIX}.x_rRNA.${Genome}.bed12 && \
	samtools sort -@ $CPU ${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam ${PREFIX}.x_rRNA.${Genome}.sorted && \
	rm -rf  ${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam && \
	samtools index ${PREFIX}.x_rRNA.${Genome}.sorted.bam && \
	paraFile=${RANDOM}.para && \
	echo "awk '\$5==255' ${PREFIX}.x_rRNA.${Genome}.bed12 | sort -k1,1 -k2,2n > ${PREFIX}.x_rRNA.${Genome}.unique.bed12 " > $paraFile && \
	echo "awk '\$5< 255' ${PREFIX}.x_rRNA.${Genome}.bed12 | sort -k1,1 -k2,2n > ${PREFIX}.x_rRNA.${Genome}.multip.bed12 " >>$paraFile && \
	ParaFly -CPU $CPU -c $paraFile -failed_cmds ${STEP}.failed_commands && \
	touch .status.${STEP}.genome_bam_processing
STEP=$((STEP+1))

###################
#Generating bigWig#
###################
paraFile="${RANDOM}${RANDOM}.para"
[ ! -f .status.${STEP}.make_bigWig ] && \
	echo "bedtools genomecov -scale $NormScale -split -bg -strand + -i ${PREFIX}.x_rRNA.${Genome}.unique.bed12 -g $Chrome > ${PREFIX}.x_rRNA.${Genome}.sorted.unique.plus.bedgraph && bedGraphToBigWig ${PREFIX}.x_rRNA.${Genome}.sorted.unique.plus.bedgraph $Chrome ${PREFIX}.x_rRNA.${Genome}.unique.plus.BW"  >  $paraFile && \
	echo "bedtools genomecov -scale $NormScale -split -bg -strand - -i ${PREFIX}.x_rRNA.${Genome}.unique.bed12 -g $Chrome | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${PREFIX}.x_rRNA.${Genome}.sorted.unique.minus.bedgraph && bedGraphToBigWig ${PREFIX}.x_rRNA.${Genome}.sorted.unique.minus.bedgraph $Chrome ${PREFIX}.x_rRNA.${Genome}.unique.minus.BW" >> $paraFile && \
	ParaFly -c $paraFile -CPU 2 && \
	echo "track name=U.${PREFIX}(+) description=${PREFIX}(+) maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.unique.plus.BW"  >  ${PREFIX}.x_rRNA.${Genome}.unique.track && \
	echo "track name=U.${PREFIX}(-) description=${PREFIX}(-) maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.unique.minus.BW" >> ${PREFIX}.x_rRNA.${Genome}.unique.track && \
	touch .status.${STEP}.make_bigWig
STEP=$((STEP+1))

#######################
#Mapping to transposon#
#######################
# mapping to transposon directly
[ ! -f .status.${STEP}.transposon_mapping_with_BWA ] && \
bwa mem \
	-t $CPU \
	-k 19 \
	-w 100 \
	-d 100 \
	-r 1.0 \
	-c 1000 \
	$BWA_TRANSPOSON_INDEX \
	${LEFT}  \
	> ${PREFIX}.x_rRNA.rep.BWA.sam && \
	touch .status.${STEP}.transposon_mapping_with_BWA
STEP=$((STEP+1))

[ ! -f .status.${STEP}.transposon_mapping_with_BWA2 ] && \
	samtools view -bS -F0x4  ${PREFIX}.x_rRNA.rep.BWA.sam  >  ${PREFIX}.x_rRNA.rep.BWA.bam && \
	rm -rf  ${PREFIX}.x_rRNA.rep.BWA.sam  && \
	samtools sort -@ $CPU ${PREFIX}.x_rRNA.rep.BWA.bam ${PREFIX}.x_rRNA.rep.BWA.sorted && \
	rm -rf ${PREFIX}.x_rRNA.rep.BWA.bam && \
touch .status.${STEP}.transposon_mapping_with_BWA2

#########################
#Mapping to piRNAcluster#
#########################
# mapping to piRNA cluster directly
[ ! -f .status.${STEP}.piRNAcluster_mapping_with_BWA ] && \
bwa mem \
	-t $CPU \
	-k 19 \
	-w 100 \
	-d 100 \
	-r 1.0 \
	-c 1000 \
	-a \
	$BWA_PIRNACLUSTER_INDEX \
	${LEFT}  \
	> ${PREFIX}.x_rRNA.cluster.BWA.sam && \
	touch .status.${STEP}.piRNAcluster_mapping_with_BWA
STEP=$((STEP+1))

[ ! -f .status.${STEP}.piRNAcluster_mapping_with_BWA2 ] && \
	samtools view -bS -f0x2  ${PREFIX}.x_rRNA.cluster.BWA.sam  >  ${PREFIX}.x_rRNA.cluster.BWA.bam && \
	rm -rf ${PREFIX}.x_rRNA.cluster.BWA.sam && \
	samtools sort -@ $CPU ${PREFIX}.x_rRNA.cluster.BWA.bam ${PREFIX}.x_rRNA.cluster.BWA.sorted && \
	rm -rf ${PREFIX}.x_rRNA.cluster.BWA.bam && \
touch .status.${STEP}.piRNAcluster_mapping_with_BWA2




