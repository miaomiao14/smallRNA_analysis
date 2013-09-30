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

INDIR=$1
OUTDIR=$2
##correlation analysis: RSQ vs DEG; RSQ vs SRA; DEG vs SRA
#FBgn0000004_17.6 in SRA, FBgn0000004_17 in DEG and RSQ

declare -a GROUPGT=("ago3cdmut" "ago3wtmut" "ago3cdwt" "aubcdmut" "aubwtmut" "aubcdwt")
declare -a ago3cdmut=("aubvasAgo3CDrescue" "ago3MutsWW")
declare -a ago3wtmut=("aubvasAgo3WTrescue" "ago3MutsWW")
declare -a ago3cdwt=("aubvasAgo3CDrescue" "aubvasAgo3WTrescue")
declare -a aubcdmut=("AubCDrescue" "aubMutsWW")
declare -a aubwtmut=("AubWTrescue" "aubMutsWW")
declare -a aubcdwt=("AubCDrescue" "AubWTrescue")

#declare -a GT=("AubCDrescue" "AubWTrescue" "aubMutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "ago3MutsWW")


#process the smRNA transposon list
#OUTDIR1=${INDIR}/smRNA/diff_nalysis
#[ ! -d ${OUTDIR1} ] && mkdir -p ${OUTDIR1}
#for g in "${GROUPGT[@]}"
#do
#	SUBGROUP="$g[@]"
#	[ -f ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.transposonpiRNAs.list.txt ] && rm ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.transposonpiRNAs.list.txt
#	[ -f ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.S.transposonpiRNAs.list.txt ] && rm ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.S.transposonpiRNAs.list.txt
#	[ -f ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.AS.transposonpiRNAs.list.txt ] && rm ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.AS.transposonpiRNAs.list.txt
#
#	for t in ${!SUBGROUP}
#	do
#		SRA=${INDIR}/smRNA/jia_pipeline_results/Phil.SRA.${t}.ox.ovary.inserts/output/Phil.SRA.${t}.ox.ovary.inserts.transposon.list
#		
#		normFactor=`cat ${INDIR}/smRNA/jia_pipeline_results/Phil.SRA.${t}.ox.ovary.inserts/output/Phil.SRA.${t}.ox.ovary.inserts_stats_table_reads|tail -1|cut -f4`
#		sed 1d $SRA| awk -v gt=$t -v nf=$normFactor '{OFS="\t"}{print gt,$1,$10/nf*1000000}' >> ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.transposonpiRNAs.list.txt
#		sed 1d $SRA| awk -v gt=$t -v nf=$normFactor '{OFS="\t"}{print gt,$1,$11/nf*1000000}' >> ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.S.transposonpiRNAs.list.txt
#		sed 1d $SRA| awk -v gt=$t -v nf=$normFactor '{OFS="\t"}{print gt,$1,$12/nf*1000000}' >> ${OUTDIR1}/Phil.SRA.${g}.ox.nnc.normlized.AS.transposonpiRNAs.list.txt
#
#	done
#done




OUTDIR2=${INDIR}/smRNA/total_sense_sense_correlation_analysis
[ ! -d ${OUTDIR2} ] && mkdir -p ${OUTDIR2}
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/A
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/M
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/RAW
[ ! -d ${OUTDIR2}/SRA_VS_DEG/A ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/A
[ ! -d ${OUTDIR2}/SRA_VS_DEG/M ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/M
[ ! -d ${OUTDIR2}/SRA_VS_DEG/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/RAW
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/A
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/M
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/RAW
paraFile=${OUTDIR2}/RSQDEGSRA.correlationanalysis.${RANDOM}.para && \
for g in "${GROUPGT[@]}"
do
	OUTDIR3=${OUTDIR2}/${g}
	SRA=${INDIR}/smRNA/diff_nalysis/Phil.SRA.${g}.ox.nnc.normlized.transposonpiRNAs.list.txt
	RSQ=${INDIR}/rnaseq/diff_analysis2/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.out.DESeqNF.transposon.normalizedcounts.txt
	DEG=${INDIR}/degradome/diff_analysis_S_AS/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.S.out.DESeqNF.transposon.normalizedcounts.txt
	
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_correlation $SRA $RSQ $DEG ${g} ${OUTDIR2}" >>${paraFile}	
done
ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands 


OUTDIR2=${INDIR}/smRNA/sense_sense_sense_correlation_analysis
[ ! -d ${OUTDIR2} ] && mkdir -p ${OUTDIR2}
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/A
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/M
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/RAW
[ ! -d ${OUTDIR2}/SRA_VS_DEG/A ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/A
[ ! -d ${OUTDIR2}/SRA_VS_DEG/M ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/M
[ ! -d ${OUTDIR2}/SRA_VS_DEG/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/RAW
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/A
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/M
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/RAW
paraFile=${OUTDIR2}/RSQDEGSRA.correlationanalysis.${RANDOM}.para && \
for g in "${GROUPGT[@]}"
do
	OUTDIR3=${OUTDIR2}/${g}
	SRA=${INDIR}/smRNA/diff_nalysis/Phil.SRA.${g}.ox.nnc.normlized.S.transposonpiRNAs.list.txt
	RSQ=${INDIR}/rnaseq/diff_analysis2/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.out.DESeqNF.transposon.normalizedcounts.txt
	DEG=${INDIR}/degradome/diff_analysis_S_AS/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.S.out.DESeqNF.transposon.normalizedcounts.txt
	
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_correlation $SRA $RSQ $DEG ${g} ${OUTDIR2}" >>${paraFile}	
done
ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands 


OUTDIR2=${INDIR}/smRNA/antisense_sense_sense_correlation_analysis
[ ! -d ${OUTDIR2} ] && mkdir -p ${OUTDIR2}
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/A
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/M
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/RAW
[ ! -d ${OUTDIR2}/SRA_VS_DEG/A ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/A
[ ! -d ${OUTDIR2}/SRA_VS_DEG/M ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/M
[ ! -d ${OUTDIR2}/SRA_VS_DEG/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/RAW
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/A
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/M
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/RAW
paraFile=${OUTDIR2}/RSQDEGSRA.correlationanalysis.${RANDOM}.para && \
for g in "${GROUPGT[@]}"
do
	OUTDIR3=${OUTDIR2}/${g}
	SRA=${INDIR}/smRNA/diff_nalysis/Phil.SRA.${g}.ox.nnc.normlized.AS.transposonpiRNAs.list.txt
	RSQ=${INDIR}/rnaseq/diff_analysis2/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.out.DESeqNF.transposon.normalizedcounts.txt
	DEG=${INDIR}/degradome/diff_analysis_S_AS/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.S.out.DESeqNF.transposon.normalizedcounts.txt
	
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_correlation $SRA $RSQ $DEG ${g} ${OUTDIR2}" >>${paraFile}	
done
ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands 



OUTDIR2=${INDIR}/smRNA/sense_antisense_sense_correlation_analysis
[ ! -d ${OUTDIR2} ] && mkdir -p ${OUTDIR2}
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/A
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/M
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/RAW
[ ! -d ${OUTDIR2}/SRA_VS_DEG/A ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/A
[ ! -d ${OUTDIR2}/SRA_VS_DEG/M ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/M
[ ! -d ${OUTDIR2}/SRA_VS_DEG/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/RAW
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/A
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/M
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/RAW
paraFile=${OUTDIR2}/RSQDEGSRA.correlationanalysis.${RANDOM}.para && \
for g in "${GROUPGT[@]}"
do
	OUTDIR3=${OUTDIR2}/${g}
	SRA=${INDIR}/smRNA/diff_nalysis/Phil.SRA.${g}.ox.nnc.normlized.S.transposonpiRNAs.list.txt
	RSQ=${INDIR}/rnaseq/diff_analysis2/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.out.DESeqNF.transposon.normalizedcounts.txt
	DEG=${INDIR}/degradome/diff_analysis_S_AS/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.AS.out.DESeqNF.transposon.normalizedcounts.txt
	
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_correlation $SRA $RSQ $DEG ${g} ${OUTDIR2}" >>${paraFile}	
done
ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands 

OUTDIR2=${INDIR}/smRNA/antisense_antisense_sense_correlation_analysis
[ ! -d ${OUTDIR2} ] && mkdir -p ${OUTDIR2}
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/A
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/M
[ ! -d ${OUTDIR2}/DEG_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/DEG_VS_RSQ/RAW
[ ! -d ${OUTDIR2}/SRA_VS_DEG/A ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/A
[ ! -d ${OUTDIR2}/SRA_VS_DEG/M ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/M
[ ! -d ${OUTDIR2}/SRA_VS_DEG/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_DEG/RAW
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/A ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/A
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/M ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/M
[ ! -d ${OUTDIR2}/SRA_VS_RSQ/RAW ] && mkdir -p ${OUTDIR2}/SRA_VS_RSQ/RAW
paraFile=${OUTDIR2}/RSQDEGSRA.correlationanalysis.${RANDOM}.para && \
for g in "${GROUPGT[@]}"
do
	OUTDIR3=${OUTDIR2}/${g}
	SRA=${INDIR}/smRNA/diff_nalysis/Phil.SRA.${g}.ox.nnc.normlized.AS.transposonpiRNAs.list.txt
	RSQ=${INDIR}/rnaseq/diff_analysis2/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.out.DESeqNF.transposon.normalizedcounts.txt
	DEG=${INDIR}/degradome/diff_analysis_S_AS/${g}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.AS.out.DESeqNF.transposon.normalizedcounts.txt
	
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_correlation $SRA $RSQ $DEG ${g} ${OUTDIR2}" >>${paraFile}	
done
ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands