#! /bin/bash -x
# targets
declare -a TARGETS=("FLY_PIRNA_CLUSTER" \
	"FLY_TRANSPOSON_ALL" \
	"FLY_TRANSPOSON_GROUP1" \
	"FLY_TRANSPOSON_GROUP2" \
	"FLY_TRANSPOSON_GROUP3" \
	"FLY_TRANSPOSON_GROUP0" \
     "FLY_flyBase_GENE" \
     "FLY_flyBase_EXON" \
     "FLY_flyBase_INTRON" \
     "FLY_flyBase_5UTR" \
     "FLY_flyBase_3UTR" \
     "FLY_REPEATMASKER_DNA" \
     "FLY_REPEATMASKER_LINE" \
     "FLY_REPEATMASKER_LTR" );
    
STORAGE_DIR=/home/hanb/scratch/CD
WD=$PWD

for t in ${TARGETS[@]}
do
	mkdir -p $t
	cd $t
	# do stuff
	# 1. small RNA
	#for i in `find $STORAGE_DIR -name "*unique*$t*bed2.gz" | grep -v unox`
		#	do
			#		i1=${i##*/}
			#		zcat $i > ${i1%.gz}
	#done
	# 2. degradome
	#for i in `find $STORAGE_DIR -name "*unique.$t*S.bed"`
		#	do
			#i1=${i##*/}
			#awk 'BEGIN{OFS="\t"}{if (ARGIND==1) {ct[$1"_"$2"_"$6]++} else {$4=ct[$1"_"$2"_"$6]; $5=1; print $0}}' $i $i | uniq > ${i1}2
	#done
	# ago3
	for g in Ago3 Aub;
		do
			for p in WT CD;
				do
					SMRNA_S=`ls  | grep  $g | grep $p | grep x_Hairpin | grep -v AS`
					SMRNA_AS=`ls | grep $g | grep $p | grep x_Hairpin | grep  AS`
					DEG_S=`ls | grep -i $g | grep $p | grep Deg | grep -v AS`
					DEG_AS=`ls | grep -i $g | grep $p | grep Deg | grep  AS`
					[ ! -z $SMRNA_S ] && [ ! -z $DEG_AS ] && echo "ppbed2 -a $SMRNA_S -b $DEG_AS > ${t}.${g}.${p}.S-AS.ppbed" >> runPP.para
					[ ! -z $SMRNA_AS ] && [ ! -z $DEG_S ] && echo "ppbed2 -a $SMRNA_AS -b $DEG_S > ${t}.${g}.${p}.AS-S.ppbed" >> runPP.para
			done
	done
	# finish 
	cd $WD
done

