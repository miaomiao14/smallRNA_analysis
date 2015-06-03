#!/bin/bash -x
# set PATH to be aware of pipeline/bin; this bin directory needs to be searched first
export PIPELINE_DIRECTORY=/diag/home/netashawang/git/RNAseq_Pipeline
export PATH=${PIPELINE_DIRECTORY}/:$PATH
#input is xkxh.norm.bed file from Jia's pipeline, gziped
NORMBED=$1
FILE=${1##*/}
inserts=`basename $FILE .norm.bed.gz`
prefix=$2
ANNOTATION_PATH=$3
OUTDIR=$4
TISSUE=$5
LOG=${inserts}.${prefix}.log
#the annotation files are in bed format
	FOLDER=$ANNOTATION_PATH
	
	#SPU_${TISSUE}_5UTR=$FOLDER/UCSC.flyBase.3UTR.bed
	SPU_${TISSUE}_FIRSTEXON=$FOLDER/${prefix}.sorted.clipped.firstexon.bed
	SPU_${TISSUE}_LASTEXON=$FOLDER/${prefix}.sorted.clipped.lastexon.bed
	SPU_${TISSUE}_MIDEXON=$FOLDER/${prefix}.sorted.clipped.middleexon.bed
	SPU_${TISSUE}_SINGLEEXON=$FOLDER/${prefix}.sorted.clipped.singleexon.bed
	SPU_${TISSUE}_ALLEXON=$FOLDER/${prefix}.sorted.clipped.allexon.bed
	SPU_${TISSUE}_INTRON=$FOLDER/${prefix}.sorted.clipped.intron.bed
	SPU_${TISSUE}_EIS=$FOLDER/${prefix}.sorted.clipped.ExonIntronSplice.slopped.bed
	SPU_${TISSUE}_IES=$FOLDER/${prefix}.sorted.clipped.IntronExonSplice.slopped.bed
	SPU_${TISSUE}_EEJ=$FOLDER/${prefix}.sorted.clipped.ExonExonJun.slopped #bowtie index
	
	declare -a TARGETS=("SPU_${TISSUE}_FIRSTEXON" \
	"SPU_${TISSUE}_LASTEXON" \
	"SPU_${TISSUE}_MIDEXON" \
    "SPU_${TISSUE}_SINGLEEXON" \
	"SPU_${TISSUE}_ALLEXON" \
	"SPU_${TISSUE}_INTRON" \
	"SPU_${TISSUE}_EIS" \
	"SPU_${TISSUE}_IES" \
	)


#normalize to the feature length, NTM, NTA, sequencing depth?
echo "$inserts" >>$LOG

[ ! -f ${OUTDIR}/$inserts.piRNA.bed ] && echo "convert the norm.bed format of piRNAs to bed format...." >>$LOG && \
zcat $NORMBED |awk -v outfile=$LOG 'BEGIN{OFS="\t"}/track/{next};{if(length($5)>=23 && length($5)<=30 ) reads=$6/$7;totalreads+=reads; print $1,$2-1,$3,$6,$7,$4,$5}END{print totalreads >outfile }' >${OUTDIR}/$inserts.piRNA.bed && \
echo "convert the norm.bed format of piRNAs to bed format done" >>$LOG
echo `date` >> $LOG
echo "mapped reads intersect to gene structures..." >> $LOG
paraFile=${OUTDIR}/${RANDOM}.para
for t in ${TARGETS[@]}
do \
echo -ne "[ ! -f ${OUTDIR}/$inserts.${t}.mapper.temp ] && bedtools intersect -a ${OUTDIR}/$inserts.piRNA.bed -b ${!t} -f 1.0 -wb > ${OUTDIR}/$inserts.${t}.mapper.temp && " >> ${paraFile} 
echo -ne "[ ! -f ${OUTDIR}/$inserts.${t}.mapper ] && awk 'BEGIN{FS=\"\\\t\";OFS=\"\\\t\"} { if(\$6==\$13){d=\"sense\";} else {d=\"antisense\";} print \$7,\$4,\$1\":\"\$2+1\"-\"\$3\"(\"\$6\")\",d,\$11,\$11,\$5}' ${OUTDIR}/$inserts.${t}.mapper.temp > ${OUTDIR}/$inserts.${t}.mapper && " >> ${paraFile}  
echo -e "[ ! -f ${OUTDIR}/$inserts.${t}.mapper2 ] && mapper2mapper2+ ${OUTDIR}/$inserts.${t}.mapper >${OUTDIR}/$inserts.${t}.mapper2" >> ${paraFile}
done

if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	CPUN=`wc -l $paraFile |cut -f1 -d" "` && \
	ParaFly -c $paraFile -CPU $CPUN -failed_cmds $paraFile.failed_commands
fi
echo `date` >> $LOG
echo "mapped reads intersect to gene structures done!" >> $LOG
paraFile=${OUTDIR}/${RANDOM}.para

#TODO: UTR EXTRACTION


#TODO:EXON-EXON JUNCTION MAPPING
##map to the exon-exon splicing junction use all raw reads
echo `date` >> $LOG
echo "raw reads mapping to the splicing junction..." >> $LOG

${PIPELINE_DIRECTORY}/run_bowtie_spu.pl $inserts.uniq.reads $mm ${SPU_${TISSUE}_EEJ} match2_EEJ.out
FLANKING=35
awk -v F=$FLANKING '{OFS="\t"; if($5>= (F-length($1)+2) && $5 <=F ) print $0;}' $inserts.match2_EEJ.out > $inserts.match2_flank35EEJ.out && \
rm $inserts.match2_EEJ.out
echo "raw reads mapping to the splicing junction done..." >> $LOG
sort -k1,1 $inserts.match2_flank35EEJ.out |bedtools groupby -g 1 -c 1 -o count >$inserts.match2_flank35EEJ.ntm   #count how many times a sequence maps to the EEJ; output format: seq NTM

${PIPELINE_DIRECTORY}/matchall2normbed.pl $inserts.match2_flank35EEJ.out $inserts.EEJ.norm.bed
normbed2mapper /home/wangw1/pipeline_spu/common/Baylor-Repeats.spurep_18.02.uniq.sort.grouped.map $inserts.norm.bed.temp map | sed 's/chr/Scaffold/g' > $inserts.transposon.mapper
mapper2mapper2+ $inserts.transposon.mapper > $inserts.transposon.mapper2

cut -f1,2 $inserts.match2_flank35EEJ.out |uniq.lines+ 0 > $inserts.match2_flank35EEJ.out.uniq.reads
uniqmap.pl $inserts.match2_flank35EEJ.out > $inserts.uniqmap.match2_flank35EEJ.out
cut -f1,2 $inserts.uniqmap.match2_flank35EEJ.out |uniq.lines+ 0 > $inserts.uniqmap.match2_flank35EEJ.out.uniq.reads
echo `date` >> $LOG
echo "raw reads mapping to the splicing junction done!" >> $LOG


for t in ${TARGETS[@]}
do \
	echo -ne " ${PIPELINE_DIRECTORY}/total_species_reads_feature_count.pl ${OUTDIR}/$inserts.${t}.mapper2 >${OUTDIR}/$inserts.${t}.mapper2.stat && " >> ${paraFile}
	echo -e  " ${PIPELINE_DIRECTORY}/genelistnew.pl ${OUTDIR}/$inserts.${t}.mapper2  >${OUTDIR}/$inserts.${t}.mapper2.genelist " >> ${paraFile}
done
if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	CPUN=`wc -l $paraFile |cut -f1 -d" "` && \
	ParaFly -c $paraFile -CPU $CPUN -failed_cmds $paraFile.failed_commands
fi