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
##correlation analysis:  DEG vs SRA
#FBgn0000004_17.6 in SRA, FBgn0000004_17 in DEG and RSQ

declare -a GT=("AubCDrescue" "AubWTrescue" "aubMutsWW" "aubvasAgo3CDrescue" "aubvasAgo3WTrescue" "ago3MutsWW")

#process the SRA transposon list, get in-cluster transposon annotation


#process the DEG transposon list


#for i in *.ovary; do a=${i}/bedIntersectWW/${i}.PE.x_rRNA.dm3.sorted.f0x40.noS.5p.all.bed.ntm.collapse.FLY_TRANSPOSON_ALL.nta.mapper2.gz;echo -e "~/git/PE_RNA_Pipeline/transposonlist.pl $a >$a.transposon.list">>parafile; done

#for i in *.ovary; do a=${i}/bedIntersectWW/${i}.PE.x_rRNA.dm3.sorted.f0x40.noS.5p.all.bed.ntm.collapse.FLY_TRANSPOSON_ALL_IN_CLUSTER.nta.mapper2.gz;echo -e "~/git/PE_RNA_Pipeline/transposonlist.pl $a >$a.transposon.list">>parafileincluster; done

OUTDIR2=${INDIR}/smRNA/senseFraction_correlation_analysis
[ ! -d ${OUTDIR2} ] && mkdir -p ${OUTDIR2}

paraFile=${OUTDIR2}/SRADEG.sensefraction.correlationanalysis.${RANDOM}.para && \
for g in "${GT[@]}"
do
	OUTDIR3=${OUTDIR2}/${g}
	SRA=${INDIR}/smRNA/jia_pipeline_results/Phil.SRA.${g}.ox.ovary.inserts/output/Phil.SRA.${g}.ox.ovary.inserts.transposon.list
	
	DEG=${INDIR}/degradome/pipeline_output_New/Phil.DEG.${g}.ovary/bedIntersectWW/Phil.DEG.${g}.ovary.PE.x_rRNA.dm3.sorted.f0x40.noS.5p.unique.bed.ntm.collapse.FLY_TRANSPOSON_ALL.nta.mapper2.gz.transposon.list
	
	echo -e "${PIPELINE_DIRECTORY}/RRR ${PIPELINE_DIRECTORY}/R.source plot_sf_correlation $SRA $DEG ${g} ${OUTDIR2}" >>${paraFile}	
done
ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands 


