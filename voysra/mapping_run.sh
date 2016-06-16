#! /bin/bash -x

# constant
genomeVersion=$2
species=$3
project=$4
BOWTIE_INDEXES=/diag/home/netashawang/scratch/bowtieIndexes/${project}/
PATH=/diag/home/netashawang/git/piPipes/bin:$PATH
ENDOGENOUS_HP_INDEX=/diag/home/netashawang/scratch/bowtieIndexes/${genomeVersion}_hairpin/${genomeVersion}.hairpin
ENDOREF=/diag/home/netashawang/scratch/bowtieIndexes/${species}.miR.hairpinCordinate.dat.bed
MM=0
CPU=8

# function delaration
function bowtie_mapping {
	_insert=$1 # insert file
	_target_name=$2 # name of target
	_target=$3 # index of target
	_prefix=${_insert%.insert}
	bowtie \
		-r \
		-S \
		-v $MM \
		-a --best --strata \
		--norc \
		-p $CPU \
		--un ${_prefix}.x_$_target_name.insert \
		$_target \
		$_insert \
		2> ${_prefix}.${_target_name}.log | \
	samtools view -bSF 0x4 - 2>/dev/null | \
	bedtools_piPipes bamtobed -i - > ${_prefix}.${_target_name}.bed && \
	piPipes_insertBed_to_bed2 $_insert ${_prefix}.${_target_name}.bed > ${_prefix}.${_target_name}.bed2 && \
	rm -rf ${_prefix}.${_target_name}.bed
	awk '{a+=$4/$5}END{printf "%.2f\n", a}' ${_prefix}.${_target_name}.bed2 > ${_prefix}.${_target_name}.count
}

# args
input_dir=$1
cd $input_dir || ( echo "cannot enter $1" && exit 1 )
input_fq=`ls | grep fastq.gz$`
prefix=${input_dir%%_*}  ##the sample name, need to be modified here
voy_hairpin_index=$BOWTIE_INDEXES/pre-${prefix}
guide_sequence=$BOWTIE_INDEXES/`ls $BOWTIE_INDEXES | grep ${input_dir%%_*} | grep -i guide` #sequence in fasta format?
passenger_sequence=$BOWTIE_INDEXES/`ls $BOWTIE_INDEXES | grep ${input_dir%%_*} | grep -i passenger`

# convert to insert file
input_insert=${input_fq%fastq.gz*}insert
if [[ ! -f $input_insert ]]; then
	piPipes_fastq_to_insert ${input_fq} $input_insert
fi

# mapping input to endogenous
echo "mapping input to endogenous miRNA"
bowtie_mapping $input_insert ${species}.hairpin $ENDOGENOUS_HP_INDEX
depth=`cat ${input_insert%.insert}.${species}.hairpin.count`
echo "mapping to artificial miRNAs"
bowtie_mapping ${input_insert%.insert}.x_${species}.hairpin.insert $prefix $voy_hairpin_index

echo "align designed sequence"
bowtie -S -f --norc $voy_hairpin_index $guide_sequence     | samtools view -bSF 0x4 - 2>/dev/null | bedtools_piPipes bamtobed -i - >  ${prefix}.designed.bed
bowtie -S -f --norc $voy_hairpin_index $passenger_sequence | samtools view -bSF 0x4 - 2>/dev/null | bedtools_piPipes bamtobed -i - >> ${prefix}.designed.bed

echo "adjust the coordinates and normalize"
bedtools_piPipes intersect -wo -s -a ${input_insert%.insert}.x_${species}.hairpin.$prefix.bed2 -b ${prefix}.designed.bed -f 0.75 | \
awk -v depth=$depth 'BEGIN{FS=OFS="\t"; nf=1000000.0/depth}{print $11,$9-$2,$3-$10,nf * $4,$5,$6,$7}' | \
sort -k1,1 -k2,2 > ${prefix}.final.bed2 && \
ln -s ${prefix}.final.bed2 final.bed2


###processing the endogenous miRNA mapping
 
bedtools_piPipes intersect -wo -a ${input_insert%.insert}.${species}.hairpin.bed2 -b $ENDOREF -f 0.75 |awk -v depth=$depth 'BEGIN{FS=OFS="\t";nf=1000000.0/depth }{print $8"-"$11,$9-$2,$3-$10,$4*nf,$5,$6,$7}' | sort -k1,1 -k2,2 >${input_insert%.insert}.${species}.hairpin.heterogeneity.bed 
## calculate N, N-1, N+1
miR.5heterogeneity.pl ${input_insert%.insert}.${species}.hairpin.heterogeneity.bed ${input_insert%.insert}.${species}.hairpin.heterogeneity.stat
