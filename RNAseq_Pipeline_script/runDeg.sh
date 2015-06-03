#! /bin/bash -x
# 6/1/2013
# Bo Han
# degradome PE library analysis pipeline

################
# Major Config #
################
export DEG_PIPELINE_VER=1.0
#whether in debug mode
export DEBUG=1 
# pipeline address
export DEG_PIPELINE_ADD=/diag/home/pkulance/scratch/PE_Pipeline
# RSEM DIRE
RSEM_DIR=/diag/home/pkulance/softwares/trinityrnaseq_r2013-02-25/trinity-plugins/rsem/
# script search PATH
export PATH=${DEG_PIPELINE_ADD}/bin:$RSEM_DIR:${PATH}
# modified bamToBed2
bamToBed2=/diag/home/pkulance/softwares/modified_tools/bedtools-2.17.0/bin/bamToBed2
#########
# USAGE #
#########
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -l left.fq -r right.fa -s fly|mouse -c 8 -o output_directory 
This is a dual library small RNA pipeline developped in the Zamore Lab in 
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

EOF
echo -en "\e[0m"
}
while getopts "hl:r:o:s:c:" OPTION
do
	case $OPTION in
		h)
			usage && exit 1
		;;
		l)
			LEFT=$OPTARG
		;;
		r)
			RIGHT=$OPTARG
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
if [[ -z $LEFT ]] || [[ -z $RIGHT ]] || [[ -z $ORGANISM ]]
then
	usage && exit 1
fi
# if CPU is undefined, then use 8
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
# if OUTDIR is undefined, then use PWD
[ ! -z ${OUTDIR} ] || OUTDIR=$PWD

##################
# Configurations #
##################
# formating date (in log)
ISO_8601='%Y-%m-%d %H:%M:%S %Z'
# step counter
STEP=1

#############
# Variables #
#############
# set ORGANISM specific files, convert upper cases into lower case so it is easier to write the case statement later
ORGANISM=`echo ${ORGANISM} | tr '[A-Z]' '[a-z]'`
# use the longest PREFIX of $LEFT and $RIGHT as the library name
LEFTFILE=${LEFT##*/}
RIGHTFILE=${RIGHT##*/}
PREFIX=`echo -e "${LEFTFILE}\n${RIGHTFILE}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && PREFIX=${PREFIX%.*}
# if $LEFT and $RIGHT does not have any PREFIX, use the name of $LEFT
[ -z "${PREFIX}" ] && PREFIX=${LEFTFILE%.f[aq]}
# Log file
LOG=${PREFIX}.${DEG_PIPELINE_VER}.log
# checking read length and determine the best sjdbOverhang
READ_LEN=`head -2 $LEFT | awk '{getline; printf "%d", length($1)}'`
# currently supported sjdbOverhang
declare -a sjdbOverhang=("49" "89" "99")
# decide which one to use
READ_LEN=`echo -e "${READ_LEN}\n${sjdbOverhang[@]}" | awk '{rl=$1;minDiff=100;choose=0;getline;for (i=1;i<=NF;++i){diff=($i-rl<0?rl-$i:$i-rl); if (diff < minDiff) {minDiff=diff; choose=$i} }}END{print choose}'`
# testing phred_format
PHRED_SCORE=`phred_test ${LEFT}`
case ${PHRED_SCORE} in
33)
	bowtiePhredOption="--phred33-quals"
	tophatPhredOption=""
;;
64)
	bowtiePhredOption="--phred64-quals"
	tophatPhredOption="--solexa1.3-quals"
;;
*)
	bowtiePhredOption=""
;;
esac

# currently, all the indexes are built with iGenome genes.gtf provided by --sjdbGTFfile and --sjdbOverhang 49/89/99
case ${ORGANISM} in
mouse)
	Genome=mm9
;;
fly)
	Genome=dm3
	BOWTIE1INDEX=/diag/home/pkulance/cloud/index/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome
	BOWTIE2INDEX=/diag/home/pkulance/cloud/index/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome
	BWAINDEX=/diag/home/pkulance/cloud/index/Drosophila_melanogaster/UCSC/dm3/Sequence/BWAIndex/genome.fa
;;
sea_urchin|seaurchin)
#TODO: sea_urchin index name and genome index building
	Genome=spu_v3.1
;;
*)
	echo "Error: unsuppported ORGANISM";
	exit 1;
;;
esac

# other variables
Chrome=${DEG_PIPELINE_ADD}/${Genome}/ChromInfo.txt
# genome fasta
GENOME_FA=${DEG_PIPELINE_ADD}/${Genome}/STARgenomeIndex/genome.fa
# directory with STAR indexes stored
GenomeIndexDir=${DEG_PIPELINE_ADD}/${Genome}/STARgenomeIndex/$READ_LEN
# directory with bowtie index for rRNA
rRNAIndex=${DEG_PIPELINE_ADD}/${Genome}/rRNAIndex/rRNA
# GTF file
GTF_FILE=${DEG_PIPELINE_ADD}/${Genome}/genes.gtf
# begin to run pipeline
[ ! -f $LOG ] && \
echo -e "`date "+$ISO_8601"`\t>beginning running PE pipeline version $DEG_PIPELINE_VER" | tee -a $LOG || \
echo -e "`date "+$ISO_8601"`\t>>resuming running PE pipeline version $DEG_PIPELINE_VER" | tee -a $LOG 

###################
# mapping to rRNA #
###################
INPUT=$rRNAIndexDir
echo -e "`date "+$ISO_8601"`\tmapping to rRNA with STAR" | tee -a $LOG
[ ! -f .status.${STEP}.rRNA_mapping ] && \
bowtie -q ${bowtiePhredOption} -v 3 \
	--fr \
	-k 1 \
	--un ${PREFIX}.x_rRNA.un.fq \
	-S \
	-p $CPU \
	${rRNAIndex} \
	-1 $LEFT \
	-2 $RIGHT \
	1>/dev/stdout \
	2>${PREFIX}.rRNA.k1.log | \
	samtools view -bSF 0x4 - > ${PREFIX}.rRNA.k1.bam && \
	touch .status.${STEP}.rRNA_mapping
	STEP=$((STEP+1))

# mapping to genome
INPUT=${GenomeIndexDir}
# redefine LEFT and RIGHT
LEFT=${PREFIX}.x_rRNA.un_1.fq
RIGHT=${PREFIX}.x_rRNA.un_2.fq
echo -e "`date "+$ISO_8601"`\tmapping to genome with STAR" | tee -a $LOG
[ ! -f .status.${STEP}.genome_mapping ] && \
	STAR \
		--runMode alignReads \
		--genomeDir $INPUT \
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
	touch .status.${STEP}.genome_mapping
STEP=$((STEP+1))
# getting basic stats
InputReads=`grep 'Number of input reads' ${PREFIX}.x_rRNA.${Genome}.Log.final.out | awk '{print $NF}'`
UniquReads=`grep 'Uniquely mapped reads number' ${PREFIX}.x_rRNA.${Genome}.Log.final.out | awk '{print $NF}'`
MultiReads=`grep 'Number of reads mapped to multiple loci' ${PREFIX}.x_rRNA.${Genome}.Log.final.out | awk '{print $NF}'`
AllMapReads=$((UniquReads+MultiReads))
UnMapReads=$((InputReads-UniquReads-MultiReads))
# normalization factor
NormScale=`echo $AllMapReads | awk '{printf "%f",1000000.0/$1}'`

#  separate unique mapper and multiple mapper
## STAR output concordant pairs only, so the NTM of left read is always the same as the right read in the same pair
echo -e "`date "+$ISO_8601"`\tprocessing mapped result" | tee -a $LOG
[ ! -f .status.${STEP}.genome_bam_processing ] && \
	samtools view -bS ${PREFIX}.x_rRNA.${Genome}.Aligned.out.sam | tee ${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam | $bamToBed2 -bedpe -split -i - > ${PREFIX}.x_rRNA.${Genome}.bedpe && \
	samtools sort -@ $CPU ${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam ${PREFIX}.x_rRNA.${Genome}.sorted && \
	rm -rf  ${PREFIX}.x_rRNA.${Genome}.Aligned.out.bam && \
	samtools index ${PREFIX}.x_rRNA.${Genome}.sorted.bam && \
	paraFile=${RANDOM}.para && \
	echo "awk '\$8==255' ${PREFIX}.x_rRNA.${Genome}.bedpe > ${PREFIX}.x_rRNA.${Genome}.unique.bedpe " > $paraFile && \
	echo "awk '\$8< 255' ${PREFIX}.x_rRNA.${Genome}.bedpe > ${PREFIX}.x_rRNA.${Genome}.multip.bedpe " >>$paraFile && \
	echo "samtools  view -b -q 255 ${PREFIX}.x_rRNA.${Genome}.sorted.bam | tee ${PREFIX}.x_rRNA.${Genome}.sorted.unique.bam | bedtools bamtobed -bed12 -split -i - | awk 'BEGIN{OFS=\"\t\"} { if (\$4~/2\$/) {\$6=(\$6==\"+\"?\"-\":\"+\");} print \$0; }' > ${PREFIX}.x_rRNA.${Genome}.sorted.unique.bed" >> $paraFile && \
	echo "samtools2 view -b -k 254 ${PREFIX}.x_rRNA.${Genome}.sorted.bam | tee ${PREFIX}.x_rRNA.${Genome}.sorted.multip.bam | bedtools bamtobed -bed12 -split -i - | awk 'BEGIN{OFS=\"\t\"} { if (\$4~/2\$/) {\$6=(\$6==\"+\"?\"-\":\"+\");} print \$0; }' > ${PREFIX}.x_rRNA.${Genome}.sorted.multip.bed" >>$paraFile && \
	ParaFly -CPU $CPU -c $paraFile -failed_cmds ${STEP}.failed_commands && \
	touch .status.${STEP}.genome_bam_processing
STEP=$((STEP+1))

# making bigWig file for both unique and multiple
echo -e "`date "+$ISO_8601"`\tmaking bigWig file" | tee -a $LOG
[ ! -f .status.${STEP}.make_bigWig ] && \
	paraFile=${RANDOM}${RANDOM}.para && \
	echo "bedtools genomecov -scale $NormScale -split -bg -strand + -i ${PREFIX}.x_rRNA.${Genome}.sorted.unique.bed -g $Chrome > ${PREFIX}.x_rRNA.${Genome}.sorted.unique.plus.bedgraph && bedGraphToBigWig ${PREFIX}.x_rRNA.${Genome}.sorted.unique.plus.bedgraph $Chrome ${PREFIX}.x_rRNA.${Genome}.unique.plus.bigWig"  >  $paraFile && \
	echo "bedtools genomecov -scale $NormScale -split -bg -strand - -i ${PREFIX}.x_rRNA.${Genome}.sorted.unique.bed -g $Chrome | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${PREFIX}.x_rRNA.${Genome}.sorted.unique.minus.bedgraph && bedGraphToBigWig ${PREFIX}.x_rRNA.${Genome}.sorted.unique.minus.bedgraph $Chrome ${PREFIX}.x_rRNA.${Genome}.unique.minus.bigWig" >> $paraFile && \
	echo "bedtools genomecov -scale $NormScale -split -bg -strand + -i ${PREFIX}.x_rRNA.${Genome}.sorted.multip.bed -g $Chrome > ${PREFIX}.x_rRNA.${Genome}.sorted.multip.plus.bedgraph && bedGraphToBigWig ${PREFIX}.x_rRNA.${Genome}.sorted.multip.plus.bedgraph  $Chrome ${PREFIX}.x_rRNA.${Genome}.multip.plus.bigWig"  >>  $paraFile && \
	echo "bedtools genomecov -scale $NormScale -split -bg -strand - -i ${PREFIX}.x_rRNA.${Genome}.sorted.multip.bed -g $Chrome | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${PREFIX}.x_rRNA.${Genome}.sorted.multip.minus.bedgraph && bedGraphToBigWig ${PREFIX}.x_rRNA.${Genome}.sorted.multip.minus.bedgraph $Chrome ${PREFIX}.x_rRNA.${Genome}.multip.minus.bigWig" >> $paraFile && \
	echo "samtools view -uf 0x42 ${PREFIX}.x_rRNA.${Genome}.sorted.unique.bam | bedtools genomecov -5 -scale $NormScale -bg -strand + -ibam stdin -g $Chrome > ${PREFIX}.x_rRNA.${Genome}.sorted.unique.5pos.plus.bedgraph   && bedGraphToBigWig ${PREFIX}.x_rRNA.${Genome}.sorted.unique.5pos.plus.bedgraph  $Chrome ${PREFIX}.x_rRNA.${Genome}.unique._1.5pos.plus.bigWig"  >> $paraFile && \
	echo "samtools view -uf 0x42 ${PREFIX}.x_rRNA.${Genome}.sorted.unique.bam | bedtools genomecov -5 -scale $NormScale -bg -strand - -ibam stdin -g $Chrome | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${PREFIX}.x_rRNA.${Genome}.sorted.unique.5pos.minus.bedgraph  && bedGraphToBigWig ${PREFIX}.x_rRNA.${Genome}.sorted.unique.5pos.minus.bedgraph $Chrome ${PREFIX}.x_rRNA.${Genome}.unique._1.5pos.minus.bigWig" >> $paraFile && \
	ParaFly -c $paraFile -CPU $CPU && \
	echo "track name=U.${PREFIX}(+) description=${PREFIX}.unique.Watson.PE maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.unique.plus.bigWig"  > ${PREFIX}.track && \
	echo "track name=U.${PREFIX}(-) description=${PREFIX}.unique.Crick.PE  maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.unique.minus.bigWig" >> ${PREFIX}.track && \
	echo "track name=M.${PREFIX}(+) description=${PREFIX}.multip.Watson.PE maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.multip.plus.bigWig"  >> ${PREFIX}.track && \
	echo "track name=M.${PREFIX}(-) description=${PREFIX}.multip.Crick.PE  maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.multip.minus.bigWig" >> ${PREFIX}.track && \
	echo "track name=U5.${PREFIX}(+) description=${PREFIX}.unique.5Pos.Watson.PE maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.unique._1.5pos.plus.bigWig"  >> ${PREFIX}.track && \
	echo "track name=U5.${PREFIX}(-) description=${PREFIX}.unique.5Pos.Crick.PE  maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${PREFIX}.x_rRNA.${Genome}.unique._1.5pos.minus.bigWig" >> ${PREFIX}.track && \
	touch .status.${STEP}.make_bigWig
STEP=$((STEP+1))

# make intersect with transposons
echo -e "`date "+$ISO_8601"`\tdoing intersecting" | tee -a $LOG
[ ! -f .status.${STEP}.intersectAll ] && \
	intersect_all.sh \
		${PREFIX}.x_rRNA.${Genome}.unique.bedpe \
		${PREFIX}.x_rRNA.${Genome}.multip.bedpe \
		${ORGANISM} \
		$GENOME_FA && \
	touch .status.${STEP}.intersectAll
STEP=$((STEP+1))

# mapping to transposon concensus sequence directly
echo -e "`date "+$ISO_8601"`\tmapping to concensus sequence of transposon directly" | tee -a $LOG
INPUT=${DEG_PIPELINE_ADD}/${Genome}/repBase/$READ_LEN
[ ! -f .status.${STEP}.transposon_mapping ] && \
	STAR \
		--runMode alignReads \
		--genomeDir $INPUT \
		--readFilesIn ${LEFT} ${RIGHT} \
		--runThreadN $CPU \
		--outFilterScoreMin 0 \
		--outFilterScoreMinOverLread 0.72 \
		--outFilterMatchNmin 0 \
		--outFilterMatchNminOverLread 0.72 \
		--outFilterMultimapScoreRange 1 \
		--outFilterMultimapNmax -1 \
		--outFilterMismatchNmax 18 \
		--outFilterMismatchNoverLmax 0.1 \
		--alignIntronMax 1 \
		--alignIntronMin 2 \
		--genomeLoad NoSharedMemory \
		--outFileNamePrefix ${PREFIX}.x_rRNA.rep. \
		--outSAMunmapped None \
		--seedSearchStartLmax 20 \
		--seedSearchStartLmaxOverLread 1.0 \
		--chimSegmentMin 0 && \
	samtools view -bS ${PREFIX}.x_rRNA.rep.Aligned.out.sam > ${PREFIX}.x_rRNA.rep.Aligned.out.bam && \
	samtools sort -@ $CPU    ${PREFIX}.x_rRNA.rep.Aligned.out.bam ${PREFIX}.x_rRNA.rep.sorted && \
	samtools sort -@ $CPU -n ${PREFIX}.x_rRNA.rep.Aligned.out.bam ${PREFIX}.x_rRNA.rep.sortedName && \
	bedtools bamtobed -i ${PREFIX}.x_rRNA.rep.sortedName.bam -bedpe > ${PREFIX}.x_rRNA.rep.sortedName.bedpe && \
	touch .status.${STEP}.transposon_mapping
STEP=$((STEP+1))


# mapping to genome with TOPHAT
LEFT=${PREFIX}.x_rRNA.un_1.fq
RIGHT=${PREFIX}.x_rRNA.un_2.fq
[ ! -f .status.${STEP}.tophat_mapping_genome ] && \
	tophat \
		-N 5 \
		--read-edit-dist 5 \
		-o ${PREFIX}.x_rRNA.tophat \
		-r 75 \
		--mate-std-dev 100 \
		$tophatPhredOption \
		-p $CPU \
		-g 100 \
		--no-discordant \
		--no-mixed \
		--no-coverage-search \
		--library-type fr-secondstrand	\
		--b2-sensitive \
		-G $GTF_FILE \
		-x 10 \
		--no-convert-bam \
		$BOWTIE2INDEX \
		$LEFT \
		$RIGHT && \
	touch .status.${STEP}.tophat_mapping_genome
STEP=$((STEP+1))






