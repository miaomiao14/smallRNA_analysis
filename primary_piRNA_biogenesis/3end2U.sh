#! /bin/bash -x



DEG_PIPELINE_ADD=/home/ww74w/bin/
BED=$1
BED_NAME=${1##*/}
Genome=$2
P=$3
OUTDIR=/home/wangw1/data/projects/cd_ip/bylength

SMRNA_PIPELINE_DIRECTORY=/home/ww74w/git/smallRNA_analysis/

##separate the piRNAs assigned to Aub, Ago3, or PIWI (from each IP) by length
##chrX	3302399	3302424	1	1	+	TTTTTTTTTTTTTTTTTTTTTAGGA

#only use unqiely mapping piRNA here
#extract the last nucleotide of each species
#evaluate the pattern TT

for i in {23..29}
do
	awk -v l=$i 'BEGIN{OFS="\t"} {if(length($7)==l) {print $0}}' ${BED} >${OUTDIR}/${BED_NAME}.$i
	###get the downstream intervals; only use unique mapping reads first
	awk 'BEGIN{OFS="\t"} { if ($5==1) { if ($6=="+") { print $1,$3,$3+20,$4,$5,$6,$7 } else { print $1,$2-20,$2,$4,$5,$6,$7 } }}' ${OUTDIR}/${BED_NAME}.$i >${OUTDIR}/${BED_NAME}.$i.uniq.+20
	bedtools nuc -bed ${OUTDIR}/${BED_NAME}.$i.uniq.+20 -fi $Genome -s -seq -pattern TT -C |tr 'tT' 'uU' >${OUTDIR}/${BED_NAME}.$i.uniq.+20.nuc.Up
	~/git/smallRNA_analysis/primary_piRNA_biogenesis/3end2U.pl ${OUTDIR}/${BED_NAME}.$i.uniq.+20.nuc.Up $P 16
	~/git/smallRNA_analysis/primary_piRNA_biogenesis/3end2U.pl ${OUTDIR}/${BED_NAME}.$i.uniq.+20.nuc.Up UU 16
	
done
