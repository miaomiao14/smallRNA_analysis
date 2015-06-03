#! /bin/bash -x
# Paired-end RNASeq pipeline in the Zamore Lab
# 2013-06-12
# Bo W Han (bo.han@umassmed.edu)
# Phillip Zamore Lab
# RNA Therapeutics Institute
# HHMI & University of Massachusetts Medical School

################
# Major Config #
################
# VERSION CONTROL
export PE_PIPELINE_VER=1.0
# IN DEBUG MODE
export DEBUG=1 
# PIPELINE DIRECTORY
export DEG_PIPELINE_ADD=/home/hanb/nearline/PE_Pipeline/
#RUN DIRECTORY
export RUN_PIPELINE_ADD=/home/wangw1/git/PE_RNA_Pipeline
# small RNA pipeline address (some scripts used)
export SMALL_RNA_PIPELINE_ADD=/home/hanb/nearline/small_RNA_Pipeline
# SCRIPTS SEARCHING PATH
export PATH=${DEG_PIPELINE_ADD}/bin:${RUN_PIPELINE_ADD}:${PATH}
# modified bamToBed2, which keep \1 and \2 order when converting bam to bedpe. a.k.a, \1 is always at col 1-3 and \2 is always 4-6
# the reason not using samtools sort + original bedtools bamtobed is because samtools sort won't perform as expected when there are secondary aligments
# and orginal bedtools bamtobed will mistakenly put adjacent bam entries, which could be primary and secondary alignments of the same end, into one bedpe.
# this modified version directly use sam output of STAR, which put proper alignments adjacent to eath other. 
BAMTOBEDPE=/home/hanb/software/modified_bedtools/bedtools-2.17.0/bin/bamToBed 

#########
# USAGE #
#########
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -l left.fq -r right.fa -s fly|mouse -c 8 -o output_directory 
This is a paired end RNA Sequencing pipeline developped in the Zamore Lab in 
University of Massachusetts Medical School. 
Please email Bo.Han@umassmed.edu for any questions or bugs. 
Thanks for using it. 

OPTIONS:
	-h      Show this message
	-l      Path to the left fastq
	-r      Path to the right fastq
	-s      Species [ fly | mouse ]
	-o      Output Directory [ default = current directory ]
	-c      Number of CPUs to use [ default = 8 ]
	-F      Library was constructed by ligation method (\1 reads is in the same direction as the transcript)
	-R      Library was constructed by dUTR method (\2 read is in the same direction as the transcript), this is the option by default
	-f      Only use the five prime end for making bigWig file and intersecting check, like CAGE and degradome; this option is false by default
EOF
echo -en "\e[0m"
}
# initiate strandless
STRAND=0
CAGE_LIKE=0
while getopts "hfl:r:o:s:c:FR" OPTION
do
	case $OPTION in
		h)
			usage && exit 1
		;;
		l)
			LEFT=`readlink -f $OPTARG`
		;;
		r)
			RIGHT=`readlink -f $OPTARG`
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
		F)
			STRAND=$((STRAND+2))
		;;
		R)
			STRAND=$((STRAND+1))
		;;
		f)
			CAGE_LIKE=1
		;;
		?)
			usage && exit 1
		;;
	esac
done
# if LEFT or RIGHT or ORGANISM is undefined, print out help and exit
if [[ -z $LEFT ]] || [[ -z $RIGHT ]] || [[ -z $ORGANISM ]]
then
	usage && exit 1
fi
# if CPU is undefined, then use 8
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
# check input
[ ! -s $LEFT ]  && echo -e "\e[1;31mError: cannot open $LEFT \e[0m" && exit 1
[ ! -s $RIGHT ] && echo -e "\e[1;31mError: cannot open $RIGHT\e[0m" && exit 1
# OUTPUT directory defined? otherwise use current directory
[ ! -z $OUTDIR ] || OUTDIR=$PWD
# test wether we need to create the new directory
[ ! -z $OUTDIR ] && mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the direcory of input fastq file\e[0m"
# enter destination direcotry, since later it will be ran in this directory, there is NO need to add OUTPUT DIR in the following codes
cd ${OUTDIR} || (echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" && exit 1)
# test writtability
touch .writting_permission && rm -rf .writting_permission || (echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" && exit 1)

# if unspecified strandless, then use dUTR method
if [[ $STRAND == 3 ]]
then
	echo -e "\e[1;31mError: cannot specify both -F and -R\e[0m"
	usage && exit 2
fi
# if not specified, then use dUTR
if [[ $STRAND == 0 ]]
then
	STRAND=1
fi

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
checkExist "htseq-count" 
checkExist "bowtie2" # mapping tool, used for most mapping purpose
checkExist "cufflinks" # mapping tool, used for most mapping purpose
checkExist "ParaFly" # a handy tool in trinity package to wrap a serials independent commands and run them with multiple threads
checkExist "intersectPE.sh" # using bedtools bedtopair to assign eash cordant pair to genomic features
checkExist "bedPE2weblogo.sh" # making information content graph
checkExist "bedSE2weblogo.sh" # making information content graph
checkExist "quantify_by_cuff.sh"
checkExist "Calculate_fpkm.sh"
checkExist "htseqCountPE.sh"
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
LEFTFILE=${LEFT##*/}
RIGHTFILE=${RIGHT##*/}
PREFIX=`echo -e "${LEFTFILE}\n${RIGHTFILE}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && PREFIX=${PREFIX%.*}
# if $LEFT and $RIGHT does not have any PREFIX, use the name of $LEFT
[ -z "${PREFIX}" ] && PREFIX=${LEFTFILE%.f[aq]}
# Log file
LOG=${OUTDIR}/${PREFIX}.${PE_PIPELINE_VER}.log
# checking read length and determine the best sjdbOverhang
READ_LEN=`head -2 $LEFT | awk '{getline; printf "%d", length($1)}'`
# overhang
STARsjdbOverhang=$((READ_LEN-1))
# testing phred_format
PHRED_SCORE=`phred_test ${LEFT}`
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
	GENOME_FA=/home/hanb/nearline/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa
;;
fly)
	Genome=dm3
	GENOME_FA=/home/hanb/nearline/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa
;;
*)
	echo -e "\e[1;31mError: Currently only mm9 and dm3 are supported\e[0m"
	usage && exit 1
;;
esac

# other variables
Chrome=${DEG_PIPELINE_ADD}/common_files/${Genome}/ChromInfo.txt
# genome fasta
GENOME_FA=${DEG_PIPELINE_ADD}/common_files/${Genome}/STARgenomeIndex/genome.fa
# directory with STAR indexes stored
STAR_Genome_INDEX_DIR=${DEG_PIPELINE_ADD}/common_files/${Genome}/STARgenomeIndex
# directory with bowtie index for rRNA
STAR_rRNA_INDEX_DIR=${DEG_PIPELINE_ADD}/common_files/${Genome}/rRNAIndex
# GTF file
GTF_FILE=${DEG_PIPELINE_ADD}/common_files/${Genome}/genes.gtf
GTF_FILE_WITH_TRANSPOSON=${DEG_PIPELINE_ADD}/common_files/${Genome}/genes_with_transposon.gff
# BWA index for transposon sequences
BWA_TRANSPOSON_INDEX=${DEG_PIPELINE_ADD}/common_files/${Genome}/repBase/transposon.ref
# BWA index for piRNA cluster
BWA_PIRNACLUSTER_INDEX=${DEG_PIPELINE_ADD}/common_files/${Genome}/piRNAcluster/piRNAcluster.fa

#######################
#Begin to run pipeline#
#######################
[ ! -f $LOG ] && \
echo -e "`date "+$ISO_8601"`\tbeginning running PE pipeline version $PE_PIPELINE_VER"  | tee -a $LOG || \
echo -e "`date "+$ISO_8601"`\t..resuming running PE pipeline version $PE_PIPELINE_VER" | tee -a $LOG 

# test if STAR index exist, if not, creat it
STAR_Genome_INDEX=$STAR_Genome_INDEX_DIR/$STARsjdbOverhang
[ -d $STAR_Genome_INDEX ] || \
	mkdir -p $STAR_Genome_INDEX || \
	STAR --runMode genomeGenerate \
		--runThreadN $CPU \
		--genomeDir $STAR_Genome_INDEX \
		--genomeFastaFiles $GENOME_FA \
		--sjdbOverhang ${STARsjdbOverhang}

########################
#ribosomal rRNA mapping#
########################
OUTDIR1=${OUTDIR}/xrRNA
[ ! -z ${OUTDIR1} ] && mkdir -p ${OUTDIR1}
cd ${OUTDIR1}
echo -e "`date "+$ISO_8601"`\tmapping to rRNA with Bowtie2" | tee -a $LOG
[ ! -f ${OUTDIR}/.status.${STEP}.rRNA_mapping ] && \
bowtie2 \
	-x $STAR_rRNA_INDEX_DIR/rRNA \
	-1 $LEFT \
	-2 $RIGHT \
	-q \
	$bowtie2PhredOption \
	--very-fast \
	-k 1 \
	--no-mixed \
	--no-discordant \
	--un-conc ${OUTDIR1}/${PREFIX}.x_rRNA.UNmapped.fq \
	-p $CPU \
	-S ${OUTDIR1}/${PREFIX}.rRNA.sam \
	2> ${OUTDIR1}/${PREFIX}.rRNA.log && \
rm -rf ${OUTDIR1}/${PREFIX}.rRNA.sam && \
touch ${OUTDIR}/.status.${STEP}.rRNA_mapping
STEP=$((STEP+1))
cd $OUTDIR
################
#genome mapping#
################
# set LEFT and RIGHT to be the unmappable reads for rRNA
OUTDIR2=${OUTDIR}/genomeMapping
[ ! -z ${OUTDIR2} ] && mkdir -p ${OUTDIR2}
cd ${OUTDIR2}
LEFT=${OUTDIR1}/${PREFIX}.x_rRNA.UNmapped.1.fq
RIGHT=${OUTDIR1}/${PREFIX}.x_rRNA.UNmapped.2.fq
[ ! -f ${OUTDIR}/.status.${STEP}.genome_mapping ] && \
STAR \
	--runMode alignReads \
	--genomeDir $STAR_Genome_INDEX \
	--readFilesIn ${LEFT} ${RIGHT} \
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
touch ${OUTDIR}/.status.${STEP}.genome_mapping
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
# !! don't delete SAM file here; it will be used later for htseq-ct
echo -e "`date "+$ISO_8601"`\tprocessing mapped result" | tee -a $LOG
[ ! -f ${OUTDIR}/.status.${STEP}.genome_bam_processing ] && \
	samtools view -bS ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.Aligned.out.sam | tee ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam | $BAMTOBEDPE -bedpe -i - > ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.bedpe && \
	samtools sort -@ $CPU ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted && \
	samtools index ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.bam && \
	samtools view -hbF0x100 ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.bam > ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.F0x100.bam && \
	awk 'BEGIN{OFS="\t";FS="\t"}{if ($8==255) {$8=1; print $0 > "/dev/stdout"} else {print $0 > "/dev/stderr"}}' ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.bedpe \
	1> ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.unique.bedpe \
	2> ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.multip.bedpe.temp && \
	awk 'BEGIN{FS="\t"; OFS="\t"} {if (ARGIND==1) {ct[$7]++} else {$8=1.0/ct[$7]; print $0}}' ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.multip.bedpe.temp ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.multip.bedpe.temp > ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.multip.bedpe && \
	rm ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.multip.bedpe.temp && \
	touch ${OUTDIR}/.status.${STEP}.genome_bam_processing
STEP=$((STEP+1))

cd $OUTDIR
################
# intersecting #
################
case $CAGE_LIKE in
0) # RNASeq/ChIPSeq
	OUTDIR3=${OUTDIR}/bigWigPE
	[ ! -z ${OUTDIR3} ] && mkdir -p ${OUTDIR3}
	echo -e "`date "+$ISO_8601"`\tprepare pair end bed for bigWig" | tee -a $LOG
	[ ! -f ${OUTDIR}/.status.${STEP}.prepare_PEbed_for_bigWig ] && \
		paraFile=${OUTDIR3}/${RANDOM}.preparePEbedforBigWig.para && \
		echo "samtools  view -b -q 255 ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.bam | tee ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.unique.bam | bedtools bamtobed -bed12 -split -i - | awk 'BEGIN{OFS=\"\t\"} { if (\$4~/${STRAND}\$/) {\$6=(\$6==\"+\"?\"-\":\"+\");} print \$0; }' > ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.unique.bed" >> $paraFile && \
		echo "samtools2 view -b -k 254 ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.bam | tee ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.multip.bam | bedtools bamtobed -bed12 -split -i - | awk 'BEGIN{OFS=\"\t\"} { if (\$4~/${STRAND}\$/) {\$6=(\$6==\"+\"?\"-\":\"+\");} print \$0; }' > ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.multip.bed" >> $paraFile && \
		ParaFly -CPU $CPU -c $paraFile -failed_cmds ${paraFile}.failed_commands && \
		touch ${OUTDIR}/.status.${STEP}.prepare_PEbed_for_bigWig
	STEP=$((STEP+1))
	
	echo -e "`date "+$ISO_8601"`\tgenerating pair end bigWig" | tee -a $LOG
	[ ! -f ${OUTDIR}/.status.${STEP}.make_PE_bigWig ] && \
		paraFile=${OUTDIR3}/${RANDOM}.makeBigWigPE.para && \
		echo "bedtools genomecov -scale $NormScale -split -bg -strand + -i ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.unique.bed -g $Chrome > ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.unique.plus.bedgraph && bedGraphToBigWig ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.unique.plus.bedgraph $Chrome ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.unique.plus.BW"  >  $paraFile && \
		echo "bedtools genomecov -scale $NormScale -split -bg -strand - -i ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.unique.bed -g $Chrome | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.unique.minus.bedgraph && bedGraphToBigWig ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.unique.minus.bedgraph $Chrome ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.unique.minus.BW" >> $paraFile && \
		ParaFly -CPU $CPU -c $paraFile -failed_cmds ${paraFile}.failed_commands && \
		echo "track name=U.${PREFIX}(+) description=${PREFIX}(+) maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.unique.plus.BW"  >  ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.unique.track && \
		echo "track name=U.${PREFIX}(-) description=${PREFIX}(-) maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.unique.minus.BW" >> ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.unique.track && \
		touch ${OUTDIR}/.status.${STEP}.make_PE_bigWig
	STEP=$((STEP+1))
	
	OUTDIR4=${OUTDIR}/bedIntersectPE
	[ ! -z ${OUTDIR4} ] && mkdir -p ${OUTDIR4}
	OUTDIR5=${OUTDIR}/webLogoPE
	[ ! -z ${OUTDIR5} ] && mkdir -p ${OUTDIR5}
	echo -e "`date "+$ISO_8601"`\tdoing intersecting for PE" | tee -a $LOG
	[ ! -f ${OUTDIR}/.status.${STEP}.intersect_PE ] && \
		intersectPE.sh \
			${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.unique.bedpe \
			${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.multip.bedpe \
			${ORGANISM} \
			$GENOME_FA \
			$STRAND \
			$OUTDIR4 \
			$OUTDIR5 && \
		touch ${OUTDIR}/.status.${STEP}.intersect_PE
	STEP=$((STEP+1))

;;
1) # CAGE/Degradome
	# for CAGE/Degradome, we get rid of pairs whose left end has softclip
	echo -e "`date "+$ISO_8601"`\tprepare bed for bigWig" | tee -a $LOG
	[ ! -f ${OUTDIR}/.status.${STEP}.prepare_bed_for_bigWig ] && \
		samtools view -hf0x40 ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.bam | \
		filterSoftClippedReads -5 | \
		samtools view -bS - | \
		bedtools bamtobed -i - | \
		awk 'BEGIN{OFS="\t";FS="\t"}{if ($6=="+") {$3=$2+1} else {$2=$3-1} print $0}' | \
			uniq \
			>  ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.bed && \
		awk 'BEGIN{OFS="\t"}{if ($5==255) {print $0 > "/dev/stdout"} else {print $0 > "/dev/stderr"}}' ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.bed \
			1> ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.bed \
			2> ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.multip.bed && \
		touch ${OUTDIR}/.status.${STEP}.prepare_bed_for_bigWig
	STEP=$((STEP+1))

#if filter those softclipped reads, the normalization factor is different from the stat from STAR final log
		
	OUTDIR3=${OUTDIR}/bigWig
	[ ! -z ${OUTDIR3} ] && mkdir -p ${OUTDIR3}
	echo -e "`date "+$ISO_8601"`\tgenerating bigWig" | tee -a $LOG
	[ ! -f ${OUTDIR}/.status.${STEP}.make_bigWig ] && \
		paraFile=${OUTDIR3}/"${RANDOM}${RANDOM}.para" && \
		echo "bedtools genomecov -scale $NormScale -bg -strand + -i ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.bed -g $Chrome > ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.plus.bedgraph && bedGraphToBigWig ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.plus.bedgraph $Chrome ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.plus.BW"  >  $paraFile && \
		echo "bedtools genomecov -scale $NormScale -bg -strand - -i ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.bed -g $Chrome | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.minus.bedgraph && bedGraphToBigWig ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.minus.bedgraph $Chrome ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.minus.BW" >> $paraFile && \
		ParaFly -CPU $CPU -c $paraFile -failed_cmds ${paraFile}.failed_commands && \
		echo "track name=U.5p.${PREFIX}(+) description=${PREFIX}(+) maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.plus.BW"  >  ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.unique.track && \
		echo "track name=U.5p.${PREFIX}(-) description=${PREFIX}(-) maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.minus.BW" >> ${OUTDIR3}/${PREFIX}.x_rRNA.${Genome}.unique.track && \
		touch ${OUTDIR}/.status.${STEP}.make_bigWig
	STEP=$((STEP+1))
	
	OUTDIR4=${OUTDIR}/bedIntersect
	[ ! -z ${OUTDIR4} ] && mkdir -p ${OUTDIR4}
	OUTDIR5=${OUTDIR}/webLogo
	[ ! -z ${OUTDIR5} ] && mkdir -p ${OUTDIR5}
	echo -e "`date "+$ISO_8601"`\tdoing intersecting" | tee -a $LOG && \
	[ ! -f ${OUTDIR}/.status.${STEP}.intersect_5END ] && \
		${RUN_PIPELINE_ADD}/intersectSE.sh \
		${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.bed \
		${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.multip.bed \
		${OUTDIR4}/${PREFIX}.intersect.summary \
		$ORGANISM \
		$GENOME_FA \
		$CPU \
		$PREFIX \
		$OUTDIR4 \
		${OUTDIR5} && \
	touch ${OUTDIR}/.status.${STEP}.intersect_5END
	STEP=$((STEP+1))
		
######################
#build sequence index#
######################
	OUTDIR6=${OUTDIR}/degIndex
	[ ! -z ${OUTDIR6} ] && mkdir -p ${OUTDIR6}
	echo -e "`date "+$ISO_8601"`\tbuilding index" | tee -a $LOG
	[ ! -f ${OUTDIR}/.status.${STEP}.build_index1 ] && \
		paraFile=${OUTDIR6}/${RANDOM}.para && \
		echo "samtools bam2fq ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.F0x100.bam | awk '{a=substr(\$1,length(\$1)); head=\$1;getline;seq=\$1; getline;getline; qual=\$1; if (a==1) {printf \"%s\\n%s\\n+\\n%s\\n\",head,seq,qual}}' > ${OUTDIR6}/${PREFIX}.x_rRNA.${Genome}.F0x100.r1.fq" >  $paraFile && \
		echo "samtools bam2fq ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.F0x100.bam | awk '{a=substr(\$1,length(\$1)); head=\$1;getline;seq=\$1; getline;getline; qual=\$1; if (a==2) {printf \"%s\\n%s\\n+\\n%s\\n\",head,seq,qual}}' > ${OUTDIR6}/${PREFIX}.x_rRNA.${Genome}.F0x100.r2.fq" >> $paraFile && \
		ParaFly -c $paraFile -CPU $CPU && \
		touch ${OUTDIR}/.status.${STEP}.build_index1
	STEP=$((STEP+1))

	[ ! -f ${OUTDIR}/.status.${STEP}.build_index2 ] && \
		paraFile=${OUTDIR6}/${RANDOM}.para && \
		echo "catFastqUnique ${OUTDIR6}/${PREFIX}.x_rRNA.${Genome}.F0x100.r1.fq > ${OUTDIR6}/${PREFIX}.x_rRNA.${Genome}.F0x100.r1.fa && bowtie-build ${OUTDIR6}/${PREFIX}.x_rRNA.${Genome}.F0x100.r1.fa ${OUTDIR6}/${PREFIX}.r1 " >  $paraFile && \
		echo "catFastqUnique ${OUTDIR6}/${PREFIX}.x_rRNA.${Genome}.F0x100.r2.fq > ${OUTDIR6}/${PREFIX}.x_rRNA.${Genome}.F0x100.r2.fa && bowtie-build ${OUTDIR6}/${PREFIX}.x_rRNA.${Genome}.F0x100.r2.fa ${OUTDIR6}/${PREFIX}.r2 " >> $paraFile && \
		ParaFly -c $paraFile -CPU $CPU && \
		echo ${OUTDIR6}/${PREFIX}.r1 >  ${OUTDIR}/indexName && \
		echo ${OUTDIR6}/${PREFIX}.r2 >> ${OUTDIR}/indexName && \
		touch ${OUTDIR}/.status.${STEP}.build_index2
	STEP=$((STEP+1))

##WW's intersect
		
	OUTDIR7=${OUTDIR}/bedIntersectWW
	[ ! -z ${OUTDIR7} ] && mkdir -p ${OUTDIR7}
	echo -e "`date "+$ISO_8601"`\tdoing intersecting" | tee -a $LOG && \
	[ ! -f ${OUTDIR}/.status.${STEP}.intersect_5END_WW ] && \
		${RUN_PIPELINE_ADD}/intersectSE_WW.sh \
		${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.unique.bed \
		${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.f0x40.noS.5p.multip.bed \
		${OUTDIR7}/${PREFIX}.intersect.summary \
		$ORGANISM \
		$GENOME_FA \
		$CPU \
		$PREFIX \
		$OUTDIR7 \
	touch ${OUTDIR}/.status.${STEP}.intersect_5END_WW
	STEP=$((STEP+1))
;;
?)
	echo -e "Error: unrecognized CAGELIKE option"; exit 1;
;;
esac

#######################
#Mapping to transposon#
#######################

case $CAGE_LIKE in
0) # RNASeq/ChIPSeq
# for RNASeq, sense is the \2, thus 0x80
	SENSE_FLAG=0x82
;;
1) # CAGE/Degradome
	SENSE_FLAG=0x42
;;
esac

# mapping to transposon directly
OUTDIR7=${OUTDIR}/transposonMapping
[ ! -z ${OUTDIR7} ] && mkdir -p ${OUTDIR7}
echo -e "`date "+$ISO_8601"`\tmapping to transposon with BWA" | tee -a $LOG
[ ! -f ${OUTDIR}/.status.${STEP}.transposon_mapping_with_BWA ] && \
bwa mem \
	-t $CPU \
	-k 19 \
	-w 100 \
	-d 100 \
	-r 1.0 \
	-c 1000 \
	$BWA_TRANSPOSON_INDEX \
	${LEFT} ${RIGHT}  \
	> ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sam && \
	touch ${OUTDIR}/.status.${STEP}.transposon_mapping_with_BWA
STEP=$((STEP+1))

[ ! -f ${OUTDIR}/.status.${STEP}.transposon_mapping_with_BWA2 ] && \
	samtools view -bS -f0x2  ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sam  >  ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.bam && \
	rm -rf  ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sam  && \
	samtools sort -@ $CPU ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.bam ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sorted && \
	rm -rf ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.bam && \
	samtools view -bf$SENSE_FLAG -F0x10 ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sorted.bam \
		| bedtools bamtobed -i - -tag AS \
		> ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sorted.f0x82F0x10.bed && \
	awk -v nf=$NormScale '{if (ARGIND==1) {a[$4]++} else {b[$1]+=1.0/a[$4]}}END{for (c in b) {print c"\t"b[c]*nf}}'  \
		${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sorted.f0x82F0x10.bed \
		${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sorted.f0x82F0x10.bed \
		> ${OUTDIR7}/${PREFIX}.x_rRNA.rep.BWA.sorted.f0x82F0x10.ct && \
	touch ${OUTDIR}/.status.${STEP}.transposon_mapping_with_BWA2
STEP=$((STEP+1))

#########################
#Mapping to piRNAcluster#
#########################

# mapping to piRNA cluster directly
OUTDIR8=${OUTDIR}/BreneckeClusterMapping
[ ! -z ${OUTDIR8} ] && mkdir -p ${OUTDIR8}

echo -e "`date "+$ISO_8601"`\tmappint to piRNA clusters" | tee -a $LOG
[ ! -f ${OUTDIR}/.status.${STEP}.piRNAcluster_mapping_with_BWA ] && \
bwa mem \
	-t $CPU \
	-k 19 \
	-w 100 \
	-d 100 \
	-r 1.0 \
	-c 1000 \
	-a \
	$BWA_PIRNACLUSTER_INDEX \
	${LEFT} ${RIGHT}  \
	> ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sam && \
	touch ${OUTDIR}/.status.${STEP}.piRNAcluster_mapping_with_BWA
STEP=$((STEP+1))

[ ! -f ${OUTDIR}/.status.${STEP}.piRNAcluster_mapping_with_BWA2 ] && \
	samtools view -bS -f0x2  ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sam  >  ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.bam && \
	rm -rf ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sam && \
	samtools sort -@ $CPU ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.bam ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sorted && \
	rm -rf ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.bam && \
	samtools view -bf$SENSE_FLAG -F0x10 ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sorted.bam \
		| bedtools bamtobed -i - -tag AS \
		> ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sorted.f0x82F0x10.bed && \
	awk -v nf=$NormScale '{if (ARGIND==1) {a[$4]++} else {b[$1]+=1.0/a[$4]}}END{for (c in b) {print c"\t"b[c]*nf}}'  \
		${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sorted.f0x82F0x10.bed \
		${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sorted.f0x82F0x10.bed \
		> ${OUTDIR8}/${PREFIX}.x_rRNA.cluster.BWA.sorted.f0x82F0x10.ct && \
	touch ${OUTDIR}/.status.${STEP}.piRNAcluster_mapping_with_BWA2
STEP=$((STEP+1))

##################################
##Counts for each feature -1 #####
##################################

#use htseq-count
OUTDIR9=${OUTDIR}/htseqCount
[ ! -z ${OUTDIR9} ] && mkdir -p ${OUTDIR9}
echo -e "`date "+$ISO_8601"`\tcount by htseq-count" | tee -a $LOG

[ ! -f ${OUTDIR}/.status.${STEP}.overlap_by_htseq-count ] && \
[ ! -f ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.Aligned.out.sam ] && samtools view -h -o ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.Aligned.out.sam ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam
${RUN_PIPELINE_ADD}/htseqCountPE.sh \
	${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.Aligned.out.sam \
	$STRAND \
	${ORGANISM} \
	$OUTDIR9 \
	$CPU && \
	#rm -rf ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.Aligned.out.sam && \
	touch ${OUTDIR}/.status.${STEP}.overlap_by_htseq-count
STEP=$((STEP+1))

##################################
##Counts for each feature -2 #####
##################################

#use raw counts, bedIntersect && groupby#########
OUTDIR10=${OUTDIR}/rawCount
[ ! -z ${OUTDIR10} ] && mkdir -p ${OUTDIR10}

echo -e "`date "+$ISO_8601"`\tcount by bed intersect" | tee -a $LOG
[ ! -f ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.all.bedpe ] && \
	cat ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.unique.bedpe ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.multip.bedpe > ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.all.bedpe 

[ ! -f ${OUTDIR}/.status.${STEP}.overlap_by_intersect ] && \
${RUN_PIPELINE_ADD}/Calculate_fpkm.sh \
	${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.all.bedpe \
	$OUTDIR10 \
	${ORGANISM} \
	$READ_LEN && \
	touch ${OUTDIR}/.status.${STEP}.overlap_by_intersect
STEP=$((STEP+1))

##################################
##Counts for each feature -3 #####
##################################

case $CAGE_LIKE in
0) # RNASeq
	LIBRARY_TYPE="fr-firststrand"
;;
1) # CAGE/Degradome
	LIBRARY_TYPE="fr-secondstrand"
;;
esac

#use cufflinks#################
OUTDIR11=${OUTDIR}/cufflinks
[ ! -z ${OUTDIR11} ] && mkdir -p ${OUTDIR11}
	echo -e "`date "+$ISO_8601"`\trpkm from cufflinks" | tee -a $LOG

[ ! -f ${OUTDIR}/.status.${STEP}.quantification_by_cuff ] && \
	cufflinks \
		-o $OUTDIR11 \
		-p $CPU \
		-G $GTF_FILE_WITH_TRANSPOSON \
		-b $GENOME_FA \
		-u \
		--library-type $LIBRARY_TYPE \
		--upper-quartile-norm \
		--total-hits-norm \
		--quiet \
		--no-update-check \
		${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.bam && \
		touch ${OUTDIR}/.status.${STEP}.quantification_by_cuff
STEP=$((STEP+1))
	
# ${RUN_PIPELINE_ADD}/quantify_by_cuff.sh ${OUTDIR2}/${PREFIX}.x_rRNA.${Genome}.sorted.bam $CAGE_LIKE ${ORGANISM} $OUTDIR11 $CPU && \
	
	
	
	
	