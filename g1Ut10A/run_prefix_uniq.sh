#!/bin/sh

#11/05/2013
#WEI WANG
#for Ping-Pong method paper, analyze the data from Yuki lab
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
INDIR=/home/wangw1/data/fly/



OUTDIR=$5
LOG=${OUTDIR}/LOG.txt

Ago3=$1
Aub=$2
gt=$3
organism=$4

zcat $Ago3 |cut -f5 |awk 'BEGIN{OFS="\t"}{ print substr($1,1,16)}' |sort -u > ${Ago3}.16prefix.uniq.species
zcat $Aub |cut -f5 |awk 'BEGIN{OFS="\t"}{ print substr($1,1,16)}' |sort -u >${Aub}.16prefix.uniq.species
	
awk 'BEGIN{OFS="\t"}{print $1,1}' ${Ago3}.16prefix.uniq.species >${Ago3}.16prefix.uniq.species2 && rm ${Ago3}.16prefix.uniq.species
awk 'BEGIN{OFS="\t"}{print $1,1}' ${Aub}.16prefix.uniq.species >${Aub}.16prefix.uniq.species2 && rm ${Aub}.16prefix.uniq.species

/home/wangw1/bin/inserts_file_methods.py -a ${Ago3}.16prefix.uniq.species2 -b ${Aub}.16prefix.uniq.species2 -u

rm ${Ago3}.16prefix.uniq.species2
rm ${Aub}.16prefix.uniq.species2

cut -f1 ${Ago3}.16prefix.uniq.species2.uniqA >${Ago3}.16prefix.uniq.species2.uniqA.weed && rm ${Ago3}.16prefix.uniq.species2.uniqA
cut -f1 ${Aub}.16prefix.uniq.species2.uniqB >${Aub}.16prefix.uniq.species2.uniqB.weed && rm ${Aub}.16prefix.uniq.species2.uniqB

zcat $Ago3 | awk 'BEGIN{OFS="\t"}{ print substr($5,1,16), $0}'  >${Ago3%.gz}.16prefix
zcat $Aub | awk 'BEGIN{OFS="\t"}{ print substr($5,1,16), $0}'  >${Aub%.gz}.16prefix

Ago3name=${Ago3/IPuniq/IPprefixuniq}
Aubname=${Aub/IPuniq/IPprefixuniq}

match.pl ${Ago3}.16prefix.uniq.species2.uniqA.weed ${Ago3%.gz}.16prefix > ${Ago3name%.gz}.16prefix && rm ${Ago3%.gz}.16prefix
match.pl ${Aub}.16prefix.uniq.species2.uniqB.weed ${Aub%.gz}.16prefix > ${Aubname%.gz}.16prefix && rm ${Aub%.gz}.16prefix




cut -f2-8 ${Ago3name%.gz}.16prefix >${Ago3name%.gz} && rm ${Ago3name%.gz}.16prefix
cut -f2-8 ${Aubname%.gz}.16prefix >${Aubname%.gz} && rm ${Aubname%.gz}.16prefix


if [ $organism = "fly" ]
then

	faFile=/home/wangw1/data/common/dmel-all-chromosome-r5.5_TAS.fasta
	indexDir=/home/wangw1/data/projects/uava/fly/bowtieIndex/
	moutDir=/home/wangw1/data/projects/uava/fly/mappingOutPut/
	queryDir=/home/wangw1/data/projects/uava/fly/querySeq/
	script=/home/wangw1/git/smallRNA_analysis/g1Ut10A/pp8_q2_ww1_zscore_sep_0515_full.pl
	
	out=/home/wangw1/data/projects/uava/fly/${gt}_prefixuniq_unox_AubIP_SRA_Ago3IP_SRA_prefix16_v15
	[ ! -f $out ] && mkdir -p $out  
	AubFile=${Aubname%.gz}
	Ago3File=${Ago3name%.gz}
	
	echo "$script -i $AubFile -j $Ago3File -o ${out} -a $faFile -b $indexDir -m $moutDir -q $queryDir -n 2 -s fly -w 16 -p 16 -d 1 -f normbed"  >> Phil.trans.parafile.pp8.prefixuniq.prefix16.${gt}
fi

if [ $organism = "bombyx" ]
then
	

	faFile=/home/wangw1/data/common/silkgenome.formatted.fa
	indexDir=/home/wangw1/data/bowtieIndex/
	moutDir=/home/wangw1/data/mappingOutPut/
	queryDir=/home/wangw1/data/querySeq/
	script=/home/wangw1/git/smallRNA_analysis/g1Ut10A/pp8_q2_ww1_zscore_sep_0515_full.pl  
#BmN4
	
	
	out=/home/wangw1/data/projects/uava/${gt}_prefixuniq_BmAgo3uniq_Siwiuniq_KNOWNTE_prefix16_v15
	[ ! -f $out ] && mkdir -p $out 
	
	AubFile=${Aubname%.gz}
	Ago3File=${Ago3name%.gz}
	
	echo "$script -i $AubFile -j $Ago3File -o ${out} -a $faFile -b $indexDir -m $moutDir -q $queryDir -n 2 -s bombyx -w 16 -p 16 -d 1 -f normbed" >>Yuki.prefixuniqbound.parafile.pp8.v15  
		
	
	
fi

if [ $organism = "mouse" ]
then
	indexDir=/home/wangw1/data/projects/uava/mouse/bowtieIndex/
	moutDir=/home/wangw1/data/projects/uava/mouse/mappingOutPut/
	queryDir=/home/wangw1/data/projects/uava/mouse/querySeq/
	faFile=/home/wangw1/data/common/mm9.formatted.fasta
	script=/home/wangw1/git/smallRNA_analysis/g1Ut10A/pp8_q2_ww1_zscore_sep_0515_full.pl  
#BmN4
	
	
	out=/home/wangw1/data/projects/uava/${gt}_prefixuniq_Mili_Miwi2_rep1_prefix16_v15
	[ ! -f $out ] && mkdir -p $out 
	
	AubFile=${Aubname%.gz}
	Ago3File=${Ago3name%.gz}
	
echo "$script -i $AubFile -j $Ago3File -o ${out} -a $faFile -b $indexDir -m $moutDir -q $queryDir -n 2 -s mouse -w 16 -p 16 -d 1 -f normbed" >>Caroll.prefixuniqbound.parafile.pp8.prefix16.v15  
	
fi

