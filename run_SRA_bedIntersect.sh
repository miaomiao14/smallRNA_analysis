#!/bin/bash -x

#input is the xkxh.norm.bed files of Jia's pipeline

	PIPELINE_DIRECTORY=/home/hanb/nearline/small_RNA_Pipeline
	FOLDER=$PIPELINE_DIRECTORY/common_files/dm3/UCSC_BEDS
	
	FLY_PIRNA_CLUSTER=$FOLDER/Brennecke.pirnaCluster.bed
	FLY_PIRNA_CLUSTER_42AB=$FOLDER/Brennecke.pirnaCluster.42AB.bed
	FLY_PIRNA_CLUSTER_FLAM=$FOLDER/Brennecke.pirnaCluster.flam.bed
	
	FLY_TRN_ALL=$FOLDER/transposon.bed2
	FLY_TRN_ALL_IN_CLUSTER=$FOLDER/transposon.inCluster.bed2
	FLY_TRN_ALL_OUT_CLUSTER=$FOLDER/transposon.outCluster.bed2

	FLY_TRN_ALL_IN_CLUSTER_IN_GENE=$FOLDER/transposon.inCluster.inGene.bed2
	FLY_TRN_ALL_IN_CLUSTER_OUT_GENE=$FOLDER/transposon.inCluster.outGene.bed2
	FLY_TRN_ALL_OUT_CLUSTER_IN_GENE=$FOLDER/transposon.outCluster.inGene.bed2
	FLY_TRN_ALL_OUT_CLUSTER_OUT_GENE=$FOLDER/transposon.outCluster.outGene.bed2
	FLY_TRN_GROUP1=$FOLDER/Zamore.group1.bed
	FLY_TRN_GROUP2=$FOLDER/Zamore.group2.bed
	FLY_TRN_GROUP3=$FOLDER/Zamore.group3.bed
	FLY_TRN_GROUP0=$FOLDER/Zamore.group0.bed

	FLY_TJ=$FOLDER/tj.bed
	FLY_HET_A=$FOLDER/HeT_A.bed
	FLY_TART=$FOLDER/TART.bed
	FLY_TAHRE=$FOLDER/TAHRE.bed
	FLY_BLOOD=$FOLDER/Blood.bed
	FLY_BURDOCK=$FOLDER/Burdock.bed
	FLY_MDG3=$FOLDER/MDG3.bed
	
	FLY_cisNATs=$FOLDER/cisNATs.bed
	FLY_STRUCTURE_LOCI=$FOLDER/structured_loci.bed
	
	FLY_TRN_ALL_IN_CLUSTER_LONG=/home/wangw1/pipeline_dm3/common/transposon.group.bed6.incluster.xflam.xcluster2
	FLY_PIRNA_GERM_CLUSTER=/home/wangw1/pipeline_dm3/common/Brennecke.pirnaCluster.x2.xflam.bed
	
	declare -a TARGETS=( \

	"FLY_PIRNA_CLUSTER" \
	"FLY_PIRNA_CLUSTER_42AB" \
	"FLY_PIRNA_CLUSTER_FLAM" \
	"FLY_TRN_ALL" \
	"FLY_TRN_ALL_IN_CLUSTER" \
	"FLY_TRN_ALL_OUT_CLUSTER" \
	"FLY_TRN_ALL_IN_CLUSTER_IN_GENE" \
	"FLY_TRN_ALL_IN_CLUSTER_OUT_GENE" \
	"FLY_TRN_ALL_OUT_CLUSTER_IN_GENE" \
	"FLY_TRN_ALL_OUT_CLUSTER_OUT_GENE" \
	"FLY_TRN_GROUP1" \
	"FLY_TRN_GROUP2" \
	"FLY_TRN_GROUP3" \
	"FLY_TRN_GROUP0" \
	"FLY_TJ" \
	"FLY_HET_A" \
	"FLY_TART" \
	"FLY_TAHRE" \
	"FLY_MDG3" \
	"FLY_BLOOD" \
	"FLY_BURDOCK" \
	"FLY_cisNATs" \
	"FLY_STRUCTURE_LOCI" )
	

for i in `ls ${INDIR}/*.inserts/*.ovary.inserts.xkxh.norm.bed.gz`
do
	OUTDIR=${i%/*}
	parafly_file=${OUTDIR}/intersect.${RANDOM}.para && \
	rm -rf $parafly_file	
	for t in ${TARGETS[@]}
	do
	#remove the first line; change the format to standard bed format; and sort the coordinates 
	echo -e " zcat $i|tail -n +2 |awk 'BEGIN{OFS="\t"}{print $1,$2-1,$$3,$6,$7,$4,$5}'| bedtools sort -i -  >${i%.norm.bed.gz}.sorted.bed2 && "" >>$parafly_file
	echo -e " bedtools intersect -a ${i%.norm.bed.gz}.sorted.bed2 -b ${!t} -f 0.999 -wo |awk 'BEGIN{FS=OFS=\"\\t\"}{ \$2+=1; if (\$6==\$18) {sense=\"sense\"} else {sense=\"antisense\"}; print \$4,\$5,\$1\":\"\$2\"-\"\$3\"(\"\$6\")\",sense,\$16,\$16,1,\$NF}' | gzip > ${INTERSECTOUTDIR}/${ALL_BED_NAME}.ntm.collapse.${t}.nta.mapper2.gz"    >> $parafly_file

	if [$t = "FLY_PIRNA_CLUSTER" ] #for cluster uniq mappers
	then
		
		
	fi

	done
	
	
done

