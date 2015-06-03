#RUN DIRECTORY
export RUN_PIPELINE_ADD=/home/wangw1/git/PE_RNA_Pipeline
Genome="dm3"
declare -a TARGETS=("FLY_PIRNA_CLUSTER" \
	"FLY_TRANSPOSON_ALL" \
	"FLY_GENE_TRANSPOSON_ALL" \
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
     "FLY_REPEATMASKER_LTR" \
     "FLY_REPEATMASKER_Satellite" \
     "FLY_REPEATMASKER_Simple_repeat" \
     "FLY_REPEATMASKER_RNA" \
     "FLY_REPEATMASKER_rRNA" \
     "FLY_REPEATMASKER_Unknown" )
declare -a GROUPGT=("ago3cdmut" "ago3wtmut" "ago3cdwt" "aubcdmut" "aubwtmut" "aubcdwt")
declare -a ago3cdmut=("aubvasAgo3CDrescue" "ago3MutsWW")
declare -a ago3wtmut=("aubvasAgo3WTrescue" "ago3MutsWW")
declare -a ago3cdwt=("aubvasAgo3CDrescue" "aubvasAgo3WTrescue")
declare -a aubcdmut=("AubCDrescue" "aubMutsWW")
declare -a aubwtmut=("AubWTrescue" "aubMutsWW")
declare -a aubcdwt=("AubCDrescue" "AubWTrescue")

declare -a str=("S" "AS")

INDIR=$1
OUT=$2
TYPE=$3
NF=$4
#For htseq-count output
#htseq-count only count sense transposon RNAs
parafly_file=${OUT}/htseqcount_df_analysis.para 
rm -rf $parafly_file

for d in ${str[@]}
do 	
	for g in "${GROUPGT[@]}"
	do
		mkdir -p ${OUT}/${g}
		OUTDIR=${OUT}/${g}
		SUBGROUP="$g[@]"
		for t in ${TARGETS[@]}
		do
			
			[ -s ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out ] && rm ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out
			[ -s ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out ] && rm ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out	
			[ -s ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out.seqDep ] && rm ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out.seqDep
			[ -s ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out.seqDep ] && rm ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out.seqDep
			for s in ${!SUBGROUP}
			do
				PREFIX=Phil.${TYPE}.${s}.ovary
				if [ $TYPE = "RSQ" ]
				then
					RSQPREFIX=Phil.${TYPE}.${s}.ovary.PE
				else
					RSQPREFIX=Phil.${TYPE}.${s}.ovary
				fi
				# getting statistics
				UniquReads=`grep 'Uniquely mapped reads number' ${INDIR}/${RSQPREFIX}/genomeMapping/${RSQPREFIX}.x_rRNA.${Genome}.Log.final.out | awk '{print $NF}'`
				MultiReads=`grep 'Number of reads mapped to multiple loci' ${INDIR}/${RSQPREFIX}/genomeMapping/${RSQPREFIX}.x_rRNA.${Genome}.Log.final.out | awk '{print $NF}'`
				AllMapReads=$((UniquReads+MultiReads))
				
			echo $s,$UniquReads >>${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out.seqDep
			echo $s,$AllMapReads >>${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out.seqDep
				if [ $NF = "TRUE" ]
				then 
					
					
					# normalization factor
					uniqNormScale=`echo $UniquReads | awk '{printf "%f",1000000.0/$1}'`
					allmNormScale=`echo $AllMapReads | awk '{printf "%f",1000000.0/$1}'`
				else
					uniqNormScale=1
					allmNormScale=1
					
				fi
				
				#if it is TRANSPOSON_ALL
				if [ $t = "FLY_TRANSPOSON_ALL" ] || [ $t = "FLY_TRANSPOSON_GROUP1" ] || [ $t = "FLY_TRANSPOSON_GROUP2" ] || [ $t = "FLY_TRANSPOSON_GROUP3" ] || [ $t = "FLY_TRANSPOSON_GROUP0" ]
				then
					echo $t
							#remove the last 5 summary lines
					#uniq	
					#head -n -5 ${INDIR}/${RSQPREFIX}/htseqCount/${RSQPREFIX}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.${d}.out | \
						#	awk '{OFS="\t"}{a=split($1,ar,"\."); print ar[1],ar[2],$2}' | \
						#bedtools groupby -i - -g 1 -c 2,3 -o collapse,sum |awk '{OFS="\t"}{b=split($2,br,",");print $1,$3  }' | \
						#awk -v nf=$uniqNormScale -v gt=$s '{OFS="\t"}{print gt,$1,$2*nf}' \
									#			>>${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out
									
					#all
					head -n -5 ${INDIR}/${RSQPREFIX}/htseqCount/${RSQPREFIX}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.${d}.out | \
						awk '{OFS="\t"}{a=split($1,ar,"\."); print ar[1],ar[2],$2}' | \
						bedtools groupby -i - -g 1 -c 2,3 -o collapse,sum |awk '{OFS="\t"}{b=split($2,br,",");print $1,$3  }' | \
						awk -v nf=$allmNormScale -v gt=$s '{OFS="\t"}{print gt,$1,$2*nf}' \
									>>${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out
							
	
				else
					#uniq
					#head -n -5 ${INDIR}/${RSQPREFIX}/htseqCount/${RSQPREFIX}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.${d}.out | \
						#	awk -v nf=$uniqNormScale -v gt=$s '{OFS="\t"}{print gt,$1,$2*nf}' >>${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out 
					#All
					head -n -5 ${INDIR}/${RSQPREFIX}/htseqCount/${RSQPREFIX}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.${d}.out | \
						awk -v nf=$allmNormScale -v gt=$s '{OFS="\t"}{print gt,$1,$2*nf}' >>${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out
				fi
			done
			#echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplot.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.uniq.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE" >> $parafly_file
			echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplot.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.${t}.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE" >> $parafly_file
			
		done
		
		#All TRANSPOSONS
		#cat ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL.uniq.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_flyBase_GENE.uniq.htseqcount.nf.${d}.out >${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out
		cat ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL.htseqcount.nf.${d}.out	${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_flyBase_GENE.htseqcount.nf.${d}.out >${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.${d}.out
		#cp ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL.uniq.htseqcount.nf.${d}.out.seqDep ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out.seqDep
		cp ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL.htseqcount.nf.${d}.out.seqDep ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.${d}.out.seqDep
		#echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplotGeneTransposon.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE" >> $parafly_file
		echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplotGeneTransposon.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE transposon" >> $parafly_file
		
		#GROUP I
		cat ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_GROUP1.htseqcount.nf.${d}.out	${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_flyBase_GENE.htseqcount.nf.${d}.out >${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_GROUP1_GENE.htseqcount.nf.${d}.out
		#cp ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL.uniq.htseqcount.nf.${d}.out.seqDep ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out.seqDep
		cp ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_GROUP1.htseqcount.nf.${d}.out.seqDep ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_GROUP1_GENE.htseqcount.nf.${d}.out.seqDep
		#echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplotGeneTransposon.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE" >> $parafly_file
		echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplotGeneTransposon.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_GROUP1_GENE.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_GROUP1_GENE.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE transposon" >> $parafly_file

		#CLUSTER

		cat ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_PIRNA_CLUSTER.htseqcount.nf.${d}.out	${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_flyBase_GENE.htseqcount.nf.${d}.out >${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_PIRNA_CLUSTER_GENE.htseqcount.nf.${d}.out
		#cp ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL.uniq.htseqcount.nf.${d}.out.seqDep ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out.seqDep
		cp ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_PIRNA_CLUSTER.htseqcount.nf.${d}.out.seqDep ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_PIRNA_CLUSTER_GENE.htseqcount.nf.${d}.out.seqDep
		#echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplotGeneTransposon.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_TRANSPOSON_ALL_GENE.uniq.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE" >> $parafly_file
		echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplotGeneTransposon.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_PIRNA_CLUSTER_GENE.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_PIRNA_CLUSTER_GENE.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE cluster" >> $parafly_file
								
		#echo "/home/wangw1/src/R-2.15.3/bin/Rscript ${RUN_PIPELINE_ADD}/MAplotGeneTransposon.r ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_GENE_TRANSPOSON_ALL.htseqcount.nf.${d}.out ${OUTDIR}/${g}.x_rRNA.dm3.Aligned.out.FLY_GENE_TRANSPOSON_ALL.htseqcount.nf.${d}.out.seqDep ${OUTDIR} $TYPE" >> $parafly_file

				
												
	done
done
if [[ ! -f ${parafly_file}.completed ]] || [[ -f $parafly_file.failed_commands ]]
then
		ParaFly -c $parafly_file -CPU 8 -failed_cmds $parafly_file.failed_commands
fi

