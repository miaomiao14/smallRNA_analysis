#!/bin/bash


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

INDIR=$1 #this is the folder store all pipeline results outmost folders
OUT=/home/wangw1/isilon_temp/smRNA/transposon_bucket
LOG=${OUT}/log

declare -a NORMFACTORTYPE=("nnc" "seqDep")

STEP=1

#transposon bucket
#declare -a GROUPGT=("pago3cdmut_ox" "pago3wtmut_ox" "pago3cdwt_ox" "pago3cdmut_unox" "pago3wtmut_unox" "pago3cdwt_unox" \
#"ago3cdmut_ox" "ago3wtmut_ox" "ago3cdwt_ox" "ago3cdmut_unox" "ago3wtmut_unox" "ago3cdwt_unox" \
#"aubcdmut_ox" "aubwtmut_ox" "aubcdwt_ox" "aubcdmut_unox" "aubwtmut_unox" "aubcdwt_unox" "aubmuthet_ox" "aubhetmut_ox" "aubKDgfpKD_unox" "gfpKDaubKD_unox" \
#"qinago3muthet_ox" "qinago3muthet_unox" \
#"AubIP_aubcdwt_ox" "AubIP_aubcdwt_unox" \
#"ago3mut_cor1_ox" "ago3mut_cor1_unox" "ago3CD_cor2_ox" "ago3CD_cor2_unox" "aubmut_cor3_ox" \
#)

declare -a GROUPGT=("pago3cdmut_ox")

declare -a pago3cdmut_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a pago3wtmut_ox=("Phil.SRA.nosAgo3WTrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a pago3cdwt_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.nosAgo3WTrescue.ox.ovary.inserts")

declare -a pago3cdmut_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a pago3wtmut_unox=("Phil.SRA.nosAgo3WTrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a pago3cdwt_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.nosAgo3WTrescue.unox.ovary.inserts")

declare -a ago3cdmut_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a ago3wtmut_ox=("Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts" "Phil.SRA.ago3MutsWW.ox.ovary.inserts")
declare -a ago3cdwt_ox=("Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.ox.ovary.inserts")

declare -a ago3cdmut_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a ago3wtmut_unox=("Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts" "Phil.SRA.ago3MutsWW.unox.ovary.inserts")
declare -a ago3cdwt_unox=("Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.aubvasAgo3WTrescue.unox.ovary.inserts")

declare -a aubcdmut_ox=("Phil.SRA.AubCDrescue.ox.ovary.inserts" "Phil.SRA.aubMutsWW.ox.ovary.inserts")
declare -a aubwtmut_ox=("Phil.SRA.AubWTrescue.ox.ovary.inserts" "Phil.SRA.aubMutsWW.ox.ovary.inserts")
declare -a aubcdwt_ox=("Phil.SRA.AubCDrescue.ox.ovary.inserts" "Phil.SRA.AubWTrescue.ox.ovary.inserts")

declare -a aubcdmut_unox=("Phil.SRA.AubCDrescue.unox.ovary.inserts" "Phil.SRA.aubMutsWW.unox.ovary.inserts")
declare -a aubwtmut_unox=("Phil.SRA.AubWTrescue.unox.ovary.inserts" "Phil.SRA.aubMutsWW.unox.ovary.inserts")
declare -a aubcdwt_unox=("Phil.SRA.AubCDrescue.unox.ovary.inserts" "Phil.SRA.AubWTrescue.unox.ovary.inserts")

declare -a aubmuthet_ox=("Phil.SRA.aubMutant.ox.ovary.inserts" "Phil.SRA.aubHets.ox.ovary.inserts")
declare -a aubhetmut_ox=("Phil.SRA.aubHets.ox.ovary.inserts" "Phil.SRA.aubMutant.ox.ovary.inserts")
declare -a aubKDgfpKD_unox=("Julius.SRA.MTD_shaub.unox.ovary.inserts" "Julius.SRA.MTD_shGFP.unox.ovary.inserts")
declare -a gfpKDaubKD_unox=("Julius.SRA.MTD_shGFP.unox.ovary.inserts" "Julius.SRA.MTD_shaub.unox.ovary.inserts")

declare -a qinago3muthet_ox=("Phil.SRA.QinAgo3Muts.ox.ovary.inserts" "Phil.SRA.QinAgo3Hets.ox.ovary.inserts")
declare -a qinago3muthet_unox=("Phil.SRA.QinAgo3Muts.unox.ovary.inserts" "Phil.SRA.QinAgo3Hets.unox.ovary.inserts")

declare -a AubIP_aubcdwt_ox=("Phil.AubIP.AubCDrescue.ox.ovary.inserts" "Phil.AubIP.AubWTrescue.ox.ovary.inserts")
declare -a AubIP_aubcdwt_unox=("Phil.AubIP.AubCDrescue.unox.ovary.inserts" "Phil.AubIP.AubWTrescue.unox.ovary.inserts")

declare -a ago3mut_cor1_ox=("Phil.SRA.ago3MutsWW.ox.ovary.inserts" "Phil.SRA.ago3MutsCJ.ox.ovary.inserts")
declare -a ago3mut_cor1_unox=("Phil.SRA.ago3MutsWW.unox.ovary.inserts" "Phil.SRA.ago3MutsCJ.unox.ovary.inserts")

declare -a ago3CD_cor2_ox=("Phil.SRA.nosAgo3CDrescue.ox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.ox.ovary.inserts")
declare -a ago3CD_cor2_unox=("Phil.SRA.nosAgo3CDrescue.unox.ovary.inserts" "Phil.SRA.aubvasAgo3CDrescue.unox.ovary.inserts")

declare -a aubmut_cor3_ox=("Phil.SRA.aubMutsWW.ox.ovary.inserts" "Phil.SRA.aubMutant.ox.ovary.inserts")

echo -e "`date` "+$ISO_8601"\tDraw paired abundance,sense_fraction of transposon piRNAs" >> $LOG
OUTDIR1=$OUT
[ ! -d $OUTDIR1 ] && mkdir -p ${OUTDIR1}


for g in "${GROUPGT[@]}"
do
	#SUBGROUP="$g[@]"
		eval "SUBGROUP=(\"\${${g}[@]}\")"  #array in bash can not be assigned directly

		[ ! -d ${OUTDIR1}/${g} ] && mkdir ${OUTDIR1}/${g}
		outputdir=${OUTDIR1}/${g}
		inputdir=$outputdir
		inputfilename1=${SUBGROUP[0]}
		inputfilename2=${SUBGROUP[1]}
		ln -s ${INDIR}/${inputfilename1}/${inputfilename1}.xkxh.transposon.mapper2.gz ${outputdir}
		ln -s ${INDIR}/${inputfilename2}/${inputfilename2}.xkxh.transposon.mapper2.gz ${outputdir}
		samplename1b=${inputfilename1#Phil.SRA.*}
		samplename2b=${inputfilename2#Phil.SRA.*}
		samplename1=${samplename1b%*.ox.ovary.inserts}
		samplename2=${samplename2b%*.ox.ovary.inserts}
		seqdepth1=`cat ${INDIR}/${inputfilename1}/output/${inputfilename1}_stats_table_reads|tail -1|awk '{print $4/1000000}'`
		seqdepth2=`cat ${INDIR}/${inputfilename2}/output/${inputfilename2}_stats_table_reads|tail -1|awk '{print $4/1000000}'`
		email="weiwanghhq@gmail.com"
		
		${PIPELINE_DIRECTORY}/bucket_new_gz_batch.pl ${inputfilename1} ${inputfilename2} ${inputdir} ${outputdir} ${samplename1} ${samplename2} ${seqdepth1} ${seqdepth2} ${email}
	
	
done


