#! /bin/bash 
# 5/2/2013
# Paired End RNA pipeline developed by Wei Wang and Bo W Han
# Zamore Lab
# University of Massachusetts Medical School and HHMI

################
# Major Config #
################
export LD_LIBRARY_PATH=$HOME/lib:$HOME/lib64
#version
export RNASEQ_PIPELINE_VER=1.0
#whether in debug mode
export DEBUG=1 

#########
# USAGE #
#########
usage() {
echo -en "\e[1;36m"
cat << EOF

usage: $0 -l left.fq -r right.fa -o output_directory -s fly|mouse -c 8 -t fr-firststrand
This pipeline is developed in the Zamore lab for the analysis of paired-end 
RNASeq made with the dUTP protocolUniversity of Massachusetts Medical School. 

Please email Wei.Wang2@umassmed.edu or Bo.Han@umassmed.edu for any questions or bugs. 
Thanks for using it. 

OPTIONS:
	-h	Show this message
	-l	Path to the left fastq
	-r	Path to the right fastq
	-s	Species [ fly | mouse ]
	-o	Output Directory [ default = current directory ]
	-c	Number of CPUs to use [ default = 8 ]
	-t	Strand type information [ fr-unstranded | fr-firststrand | fr-secondstrand ]
	-e	Pass user defined variable to the pipeline

EOF
echo -en "\e[0m"
}
while getopts "hl:r:o:s:c:t:e:" OPTION
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
	o)
		OUTDIR=$OPTARG
	;;
	s)
		ORGANISM=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'`
	;;
	c)
		CPU=$OPTARG 
	;;
	t)
		LIB_TYPE=$OPTARG
	;;
	e)
		EXTRA_VAR=${OPTARG} # using ${} to cover the situation when there is white space in OPTARG
	;;
	?)
		usage && exit 1
	;;
	esac
done
# if LEFT or RIGHT or ORGANISM is unset, exit
if [[ -z $LEFT ]] || [[ -z $RIGHT ]] || [[ -z $ORGANISM ]]
then
	usage && exit 1
fi
# if OUTDIR is not set, use current directory
OUTDIR=${OUTDIR:-$PWD}
# if CPU is not set or is not pure num
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
# if LIB_TYPE is not set
LIB_TYPE=${LIB_TYPE:-fr-unstranded}
if [ "${LIB_TYPE}" != "fr-unstranded" ] && [ "${LIB_TYPE}" != "fr-firststrand" ] && [ "${LIB_TYPE}" != "fr-secondstrand" ];
then
	echo -e "\e[1;31mError: unrecognized library type: choose from fr-unstranded | fr-firststrand | fr-secondstrand\e[0m"
	usage && exit
fi

# formating date (in log)
ISO_8601='%Y-%m-%d %H:%M:%S %Z'
# use the longest prefix of $LEFT and $RIGHT as the library name
LEFTFILE=${LEFT##*/}
RIGHTFILE=${RIGHT##*/}
prefix=`echo -e "${LEFTFILE}\n${RIGHTFILE}" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && prefix=${prefix%.*}
prefix=`echo $prefix | sed -e s/[^a-z0-9]*$//g`
# if $LEFT and $RIGHT does not have any prefix, use the name of $LEFT
[ -z "${prefix}" ] && prefix=${LEFTFILE%.f[aq]}
# Log file
LOG=${OUTDIR}/${prefix}.${RNASEQ_PIPELINE_VER}.log
# checking read length and determine the best sjdbOverhang
READ_LEN=`head -2 $LEFT | awk '{getline; printf "%d", length($1)}'`
# currently supported sjdbOverhang
declare -a sjdbOverhang=("49" "89" "99")

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
checkExist "bowtie" # mapping tool, used for most mapping purpose
checkExist "ParaFly" # a handy tool in trinity package to wrap a serials independent commands and run them with multiple threads
echo -e "\e[1;35mDone with testing required softwares/scripts, starting pipeline...\e[0m"

####################
# running pipeline #
####################
case ${ORGANISM} in
mouse)
	Genome=mm9
;;
fly)
	Genome=dm3
;;
sea_urchin|seaurchin)
	Genome=spu
;;
*)
	echo "\e[1;36mError: unsuppported organism...\e[0m" && exit 1;
;;
esac

# other variables
Chrome=${DEG_PIPELINE_ADD}/${Genome}/ChromInfo.txt
# directory with STAR indexes stored
GenomeIndexDir=${DEG_PIPELINE_ADD}/${Genome}/STARgenomeIndex/$READ_LEN
# directory with bowtie index for rRNA
rRNAIndex=${DEG_PIPELINE_ADD}/${Genome}/rRNAIndex/rRNA
# eval extra, to overwirte preset variable
eval ${EXTRA_VAR}

[ ! -f $LOG ] && \
echo -e "`date "+$ISO_8601"`\t>beginning running PE pipeline version $RNASEQ_PIPELINE_VER" >  $LOG || \
echo -e "`date "+$ISO_8601"`\t>>resuming running PE pipeline version $RNASEQ_PIPELINE_VER" >> $LOG 

#step counter
STEP=1
# TODO: make this standard
GenomeDir=/diag/home/netashawang/cloud/spu_ncRNA_star
OUTPUT1=${OUTDIR}/${prefix}_xrRNA
echo -e " step $STEP" >> $LOG
echo -e "`date "+$ISO_8601"`\tmapping to rRNA with STAR" >> $LOG 
[ ! -f ${OUTDIR}/.status.${STEP}.rRNA_mapping ] && \
mkdir $OUTPUT1 && \
cd $OUTPUT1 && \
STAR --genomeDir $GenomeDir --readFilesIn $LEFT $RIGHT --runThreadN $CPU \
	--outFilterScoreMin 0 \
	--outFilterScoreMinOverLread 0.66 \
	--outFilterMatchNmin 0 \
	--outFilterMatchNminOverLread 0.66 \
	--outFilterMultimapScoreRange 1 \
	--outFilterMultimapNmax 10000 \
	--outFilterMismatchNoverLmax 0.1 \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--outReadsUnmapped Fastx \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--genomeLoad NoSharedMemory \
	--outSAMattributes Standard \
	--limitIObufferSize 500000000 && \
echo -e "`date "+$ISO_8601"`\tmapping to rRNA with STAR done!" >> $LOG && \
echo -e "`date "+$ISO_8601"`\tchanging file names" >> $LOG && \
mv Unmapped.out.mate1 ${prefix}_xrRNA_STAR_left.fq && \
mv Unmapped.out.mate2 ${prefix}_xrRNA_STAR_right.fq && \
mv Log.final.out ${prefix}_xrRNA.Log.final.out && \
mv Log.out ${prefix}_xrRNA.Log.out && \
mv Log.progress.out ${prefix}_xrRNA.Log.progress.out && \
mv SJ.out.tab ${prefix}_xrRNA.SJ.out.tab && \
echo -e "`date "+$ISO_8601"`\tconvert sam to bam and sort..." >> $LOG && \
samtools view -bS Aligned.out.sam > Aligned.out.bam && \
samtools sort -@ 12 -m 10G Aligned.out.bam ${prefix}.STAR.sorted && \
samtools index ${prefix}.STAR.sorted.bam && \
rm Aligned.out.sam && \
rm Aligned.out.bam && \
touch ${OUTDIR}/.status.${STEP}.rRNA_mapping 
STEP=$((STEP+1))

#Mapping the leftiover reads to the genome and splicing junctions

READ_LEN=`head -2 $LEFT | awk '{getline; printf "%d", length($1)}'`
FLANKING=$((READ_LEN-1))
echo -e " step $STEP" >> $LOG
echo -e "`date "+$ISO_8601"`\tmap the non rRNA reads to genome with STAR..." >> $LOG
#Generating the star index with SJDB if it does not exist
GenomeDir_SJDB=/diag/home/netashawang/cloud/${ORGANISM}_star_SJDB_${FLANKING}
[ ! -d $GenomeDir_SJDB ] && mkdir $GenomeDir_SJDB && \
echo -e "`date "+$ISO_8601"`\tif there is no STAR index for this read length ${FLANKING}, generating a new one..." >> $LOG && \
STAR --runMode genomeGenerate --genomeDir $GenomeDir_SJDB --genomeFastaFiles /diag/home/netashawang/indexes/Spur_3.1.LinearScaffold.fa --sjdbGTFfile ${GTF} --sjdbOverhang ${FLANKING} --runThreadN $CPU --limitGenomeGenerateRAM 640000000000 && \
echo -e "`date "+$ISO_8601"`\tgenerating a new genome index with SJDB done..." >> $LOG && \
touch ${OUTDIR}/.status.${STEP}.${ORGANISM}_star_SJDB_${FLANKING}.index.generating 
STEP=$((STEP+1)) 

xrRNA_LEFT=${OUTPUT1}/${prefix}_xrRNA_STAR_left.fq
xrRNA_RIGHT=${OUTPUT1}/${prefix}_xrRNA_STAR_right.fq
OUTPUT2=$OUTDIR/${prefix}_xrRNA_genome_sjdb_mapping
cd $OUTDIR
echo -e " step $STEP" >> $LOG
echo -e "`date "+$ISO_8601"`\tmap the non rRNA reads to genome with STAR with SJDB..." >> $LOG
[ ! -f ${OUTDIR}/.status.${STEP}.genome_sjdb_mapping ] && \
mkdir $OUTPUT2 && \
cd $OUTPUT2 && \
if [[ $LIB_TYPE == "fr-unstranded" ]]
then
STAR --genomeDir $GenomeDir_SJDB --readFilesIn $xrRNA_LEFT $xrRNA_RIGHT --runThreadN $CPU \
	--outFilterScoreMin 0 \
	--outFilterScoreMinOverLread 0.66 \
	--outFilterMatchNmin 0 \
	--outFilterMatchNminOverLread 0.66 \
	--outFilterMultimapScoreRange 1 \
	--outFilterMultimapNmax 10000 \
	--outFilterMismatchNoverLmax 0.1 \
	--alignIntronMax 40000 \
	--alignIntronMin 21 \
	--alignMatesGapMax 1000000 \
	--genomeLoad NoSharedMemory \
	--outSAMattributes Standard \
	--outSAMstrandField intronMotif \
	--outReadsUnmapped Fastx \
	--limitIObufferSize 2500000000 \
	--sjdbScore 0 \
	--alignSJDBoverhangMin 8 \
	--outSJfilterCountUniqueMin 4 2 2 2 \
	--outSJfilterCountTotalMin 4 2 2 2 
else
STAR --genomeDir $GenomeDir_SJDB --readFilesIn $xrRNA_LEFT $xrRNA_RIGHT --runThreadN $CPU \
        --outFilterScoreMin 0 \
        --outFilterScoreMinOverLread 0.66 \
        --outFilterMatchNmin 0 \
        --outFilterMatchNminOverLread 0.66 \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 10000 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMax 40000 \
        --alignIntronMin 21 \
        --alignMatesGapMax 1000000 \
        --genomeLoad NoSharedMemory \
        --outSAMattributes Standard \
        --outReadsUnmapped Fastx \
        --limitIObufferSize 2500000000 \
        --sjdbScore 0 \
        --alignSJDBoverhangMin 8 \
        --outSJfilterCountUniqueMin 4 2 2 2 \
        --outSJfilterCountTotalMin 4 2 2 2 
fi
[ $? == 0 ] && \
mv Unmapped.out.mate1 ${prefix}_xrRNA_STAR_left.fq && \
mv Unmapped.out.mate2 ${prefix}_xrRNA_STAR_right.fq && \
mv Log.final.out ${prefix}_xrRNA_genome_sjdb_STARmapping.Log.final.out && \
mv Log.out ${prefix}_xrRNA_genome_sjdb_STARmapping.Log.out && \
mv Log.progress.out ${prefix}_xrRNA_genome_sjdb_STARmapping.Log.progress.out && \
mv SJ.out.tab ${prefix}_xrRNA_genome_sjdb_STARmapping.SJ.out.tab && \
samtools view -hbS Aligned.out.sam -o Aligned.out.bam && \
samtools sort -@ 12 -m 10G Aligned.out.bam ${prefix}_xrRNA_genome_sjdb_STARmapping.sorted && \
samtools index ${prefix}_xrRNA_genome_sjdb_STARmapping.sorted.bam && \
rm Aligned.out.sam && \
rm Aligned.out.bam && \
echo -e "`date "+$ISO_8601"`\tmap the non rRNA reads to genome with STAR with SJDB done!" >> $LOG && \
touch ${OUTDIR}/.status.${STEP}.genome_sjdb_mapping 
STEP=$((STEP+1))

OUTPUT3=$OUTDIR/${prefix}_cufflinks.assembly
echo -e " step $STEP " >> $LOG
echo -e " library type ${LIB_TYPE} " >> $LOG 
echo -e "`date "+$ISO_8601"`\t use cufflinks to assemble and quantitate the transcripts..." >> $LOG && \
[ ! -f ${OUTDIR}/.status.${STEP}.cufflinks.assembly ] && \
cd ${OUTDIR} && \
mkdir ${OUTPUT3} && \
cufflinks -v --no-update-check \
	--library-type ${LIB_TYPE} \
	-p ${CPU} \
	-g ${GTF} \
	-G ${GTF} \
	-o ${OUTPUT3} \
	${OUTPUT2}/${prefix}_xrRNA_genome_sjdb_STARmapping.sorted.bam \
	-u \
	-j 0.2 \
	--min-frags-per-transfrag 40 \
	--overlap-radius 100 \
	2&> ${OUTPUT3}/cufflinks.log && \
mv ${OUTPUT3}/transcripts.gtf ${OUTPUT3}/${prefix}_cufflink_transcripts.gtf && \
mv ${OUTPUT3}/genes.fpkm_tracking ${OUTPUT3}/${prefix}_cufflink_genes.fpkm_tracking && \
mv ${OUTPUT3}/isoforms.fpkm_tracking ${OUTPUT3}/${prefix}_cufflink_isoforms.fpkm_tracking && \
echo -e "`date "+$ISO_8601"`\t use cufflinks to assemble or quantitate the transcripts done" >> $LOG && \
touch ${OUTDIR}/.status.${STEP}.cufflinks.assembly
STEP=$((STEP+1))

OUTPUT4=$OUTDIR/${prefix}_cufflinks.quantification 
echo -e " step $STEP " >> $LOG 
echo -e "`date "+$ISO_8601"`\t use cufflinks to quantitate the transcripts..." >> $LOG && \
[ ! -f ${OUTDIR}/.status.${STEP}.cufflinks.quantification ] && \
cd ${OUTDIR} && \
mkdir ${OUTPUT4} && \
cufflinks -v --no-update-check \
        --library-type ${LIB_TYPE} \
        -p ${CPU} \
        -G ${GTF} \
        -o ${OUTPUT4} \
        ${OUTPUT2}/${prefix}_xrRNA_genome_sjdb_STARmapping.sorted.bam \
        -u \
        -j 0.2 \
        --min-frags-per-transfrag 40 \
	--overlap-radius 100 \
        2&> ${OUTPUT4}/cufflinks.quantification.log && \
mv ${OUTPUT4}/transcripts.gtf ${OUTPUT4}/${prefix}_cufflink_transcripts.gtf && \
mv ${OUTPUT4}/genes.fpkm_tracking ${OUTPUT4}/${prefix}_cufflink_quantification_genes.fpkm_tracking && \
mv ${OUTPUT4}/isoforms.fpkm_tracking ${OUTPUT4}/${prefix}_cufflink_quantification_isoforms.fpkm_tracking && \
echo -e "`date "+$ISO_8601"`\t use cufflinks to quantitate the transcripts done" >> $LOG && \
touch ${OUTDIR}/.status.${STEP}.cufflinks.quantification
STEP=$((STEP+1))


#  separate unique mapper and multiple mapper
## STAR output concordant pairs only, so the NTM of left read is always the same as the right read in the same pair
#TODO: need test
echo -e "`date "+$ISO_8601"`\tprocessing mapped result" >> $LOG
[ ! -f .status.${STEP}.genome_bam_processing ] && \
	samtools view -bS ${prefix}.x_rRNA.${Genome}.Aligned.out.sam > ${prefix}.x_rRNA.${Genome}.Aligned.out.bam && \
	samtools sort -@ $CPU ${prefix}.x_rRNA.${Genome}.Aligned.out.bam ${prefix}.x_rRNA.${Genome}.sorted && \
	rm -rf  ${prefix}.x_rRNA.${Genome}.Aligned.out.bam && \
	samtools index ${prefix}.x_rRNA.${Genome}.sorted.bam && \
	paraFile=${RANDOM}.para && \
	echo "samtools  view -b -q 255 ${prefix}.x_rRNA.${Genome}.sorted.bam | tee ${prefix}.x_rRNA.${Genome}.sorted.unique.bam | bedtools bamtobed -bed12 -split -i - | awk 'BEGIN{OFS=\"\t\"} { if (\$4~/1\$/) {\$6=(\$6==\"+\"?\"-\":\"+\");} print \$0; }' > ${prefix}.x_rRNA.${Genome}.sorted.unique.bed" > $paraFile && \
	echo "samtools2 view -b -k 254 ${prefix}.x_rRNA.${Genome}.sorted.bam | tee ${prefix}.x_rRNA.${Genome}.sorted.multip.bam | bedtools bamtobed -bed12 -split -i - | awk 'BEGIN{OFS=\"\t\"} { if (\$4~/1\$/) {\$6=(\$6==\"+\"?\"-\":\"+\");} print \$0; }' > ${prefix}.x_rRNA.${Genome}.sorted.multip.bed" >>$paraFile && \
	ParaFly -CPU $CPU -c $paraFile -failed_cmds ${STEP}.failed_commands && \
	touch .status.${STEP}.genome_bam_processing
STEP=$((STEP+1))

# extract flanking regions of splicing junctions
awk 'BEGIN{OFS="\t";}{if($4==1 || $4==2) print $0}' SJ.out.tab |awk 'BEGIN{OFS="\t";}  {if($4==1) {strand="+";} if($4==2){strand="-";} print $1, $2-1, $3, "STARjunc"NR, 0, strand;}'  >${prefix}.SJ.out.bed


# making bigWig file

echo -e "`date "+$ISO_8601"`\tmaking bigWig file" >> $LOG
[ ! -f .status.${STEP}.make_bigWig ] && \
	paraFile=${RANDOM}${RANDOM}.para && \
	echo "bedtools genomecov -split -bg -strand + -i ${prefix}.x_rRNA.${Genome}.sorted.unique.bed -g $Chrome > ${prefix}.x_rRNA.${Genome}.sorted.unique.plus.bedgraph  && bedGraphToBigWig ${prefix}.x_rRNA.${Genome}.sorted.unique.plus.bedgraph  $Chrome ${prefix}.x_rRNA.${Genome}.sorted.unique.plus.bigWig"  >  $paraFile && \
	echo "bedtools genomecov -split -bg -strand - -i ${prefix}.x_rRNA.${Genome}.sorted.unique.bed -g $Chrome | awk 'BEGIN{OFS=\"\t\"}{\$4 = -\$4; print \$0}' > ${prefix}.x_rRNA.${Genome}.sorted.unique.minus.bedgraph && bedGraphToBigWig ${prefix}.x_rRNA.${Genome}.sorted.unique.minus.bedgraph $Chrome ${prefix}.x_rRNA.${Genome}.sorted.unique.minus.bigWig" >> $paraFile && \
	ParaFly -c $paraFile -CPU $CPU && \
	echo "track name=${prefix}(+) description=${prefix}.unique.Watson.PE maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${prefix}.x_rRNA.${Genome}.sorted.unique.plus.bigWig"  > ${prefix}.track && \
	echo "track name=${prefix}(-) description=${prefix}.unique.Crick.PE  maxHeightPixels=25 alwaysZero=on  autoScale=on  yLineMark=0 yLineOnOff=on type=bigWig color=$((RANDOM%255)),$((RANDOM%255)),$((RANDOM%255)) visibility=full bigDataUrl=http://zlab.umassmed.edu/~hanb/Bo/${prefix}.x_rRNA.${Genome}.sorted.unique.minus.bigWig" >> ${prefix}.track && \
	touch .status.${STEP}.make_bigWig
STEP=$((STEP+1))


	
	
