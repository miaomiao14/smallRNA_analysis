#!/usr/bin/env bash

# calculate the aggregation on "metagenes" (average over all gene)
# modify from Shikui's pipeline
# 2014-08-28
# wei.wang2@umassmed.edu

infile=$1  ##norm.bed files from Jia's pipeline
in=${1##*/}
dir=${1%/$in}
inprefix=`basename $in .xkxh.norm.bed.gz`
#inprefix=${in%*.xkxh.norm.bed}


chrom_sizes=/home/wangw1/pipeline_spu/common/Spur_3.1.scaffoldsize ## 
gdb=dm3
tissue=$2
genebed=/home/tg92w/project/common/UCSC.flyBase.3UTR.bed
normFac=`grep ${inprefix} /home/wangw1/pipeline_spu/common/spu.norm.factor | cut -f2`

tssUpLen=200
ttsDnLen=200
nPtsGenebody=1000
grid=20


declare -i nPts_uptss=$tssUpLen/$grid
declare -i nPts_dntts=$ttsDnLen/$grid

zcat $infile | awk -F "\t" -v norm=$normFac '{OFS="\t"; if($4=="+" && length($5)>=23) print $1,$2-1,$3,$6","$7","norm/1000000}' > ${dir}/${inprefix}.bed.plus.pirna
zcat $infile | awk -F "\t" -v norm=$normFac '{OFS="\t"; if($4=="-" && length($5)>=23) print $1,$2-1,$3,$6","$7","norm/1000000}' > ${dir}/${inprefix}.bed.minus.pirna
## Scaffold721	25700	25729	2,1,1000

[ ! -f ${dir}/${inprefix}.bed.plus.pirna ] && \
bedClip ${dir}/${inprefix}.bed.plus.pirna $chrom_sizes ${dir}/${inprefix}.bed.plus.pirna.clip.bed && \
rm ${dir}/${inprefix}.bed.plus.pirna

[ ! -f ${dir}/${inprefix}.bed.minus.pirna ] && \
bedClip ${dir}/${inprefix}.bed.minus.pirna $chrom_sizes ${dir}/${inprefix}.bed.minus.pirna.clip.bed && \
rm ${dir}/${inprefix}.bed.minus.pirna

[ ! -f ${dir}/${inprefix}.bed.plus.pirna.clip.bed ] && \
sortBed -i ${dir}/${inprefix}.bed.plus.pirna.clip.bed | bedItemOverlapCountWithScore spu stdin -chromSize=$chrom_sizes > ${dir}/${inprefix}.plus.pirna.norm.bedGraph && \
rm ${dir}/${inprefix}.bed.plus.pirna.clip.bed

[ ! -f ${dir}/${inprefix}.bed.minus.pirna.clip.bed ] && \
sortBed -i ${dir}/${inprefix}.bed.minus.pirna.clip.bed | bedItemOverlapCountWithScore spu stdin -chromSize=$chrom_sizes > ${dir}/${inprefix}.minus.pirna.bedGraph && \
rm ${dir}/${inprefix}.bed.minus.pirna.clip.bed

[ ! -f ${dir}/${inprefix}.minus.pirna.bedGraph ] && \
awk -F "\t" '{OFS="\t"; print $1,$2,$3,-$4}' ${dir}/${inprefix}.minus.pirna.bedGraph > ${dir}/${inprefix}.minus.pirna.norm.bedGraph && \
rm ${dir}/${inprefix}.minus.pirna.bedGraph

[ ! -f ${dir}/${inprefix}.plus.pirna.norm.bedGraph ] && \
bedGraphToBigWig ${dir}/${inprefix}.plus.pirna.norm.bedGraph $chrom_sizes ${dir}/${inprefix}.plus.pirna.norm.bw && \
rm ${dir}/${inprefix}.plus.pirna.norm.bedGraph
[ ! -f ${dir}/${inprefix}.minus.pirna.norm.bedGraph ] && \
bedGraphToBigWig ${dir}/${inprefix}.minus.pirna.norm.bedGraph $chrom_sizes ${dir}/${inprefix}.minus.pirna.norm.bw && \
rm ${dir}/${inprefix}.minus.pirna.norm.bedGraph 

#Shikui's method
#sed -e '1d' $infile | awk -F "\t" {if(length($5)>=23) print $1,$2-1,$3,$6,$7,$4}
## it does not work, as it does not take the reads and NTM into consideration. The reason is he uses bam file as input

### get gene information
[ ! -f ${dir}/gene_tss_tts.bed ] && cut -f1-6 $genebed | sort -k4,4 -u | sort -k1,1 -k2,2n > ${dir}/gene_tss_tts.bed

## get gene bucket value
fetchBWscore_genebucket.pl ${dir}/${inprefix}.plus.pirna.norm.bw ${dir}/gene_tss_tts.bed $tssUpLen $nPtsGenebody $ttsDnLen $grid >${dir}/${inprefix}.plus.pirna.genebucket
fetchBWscore_genebucket.pl ${dir}/${inprefix}.minus.pirna.norm.bw ${dir}/gene_tss_tts.bed $tssUpLen $nPtsGenebody $ttsDnLen $grid >${dir}/${inprefix}.minus.pirna.genebucket

echo """ #!/home/wangw1/bin/Rscript
args <- commandArgs (TRUE);

data_file_plus <- args[1];
data_file_minus <- args[2];
outprefix <- args[3];

#read data: nGen*nWid
bkt_val_plus <- read.table (data_file_plus, sep=\"\t\", header=FALSE);
bkt_val_minus <- read.table (data_file_minus,sep=\"\t\", header=FALSE);
# average signals across rows (genes)
avg_bkt_plus <- apply ( bkt_val_plus, 2, mean);
avg_bkt_minus <- apply ( bkt_val_minus,2,mean);
# plot the avg signals
pdf( paste (outprefix, \"piRNA_genebucket.pdf\", sep=\".\"));
plot ( c((1-$nPts_uptss):(length(avg_bkt)-$nPts_uptss)), avg_bkt_plus,ylim=c(-max(avg_bkt_plus),max(avg_bkt_plus)), type=\"l\", main= \"gene.aggplot.up$tssUpLen.down$ttsDnLen.${grid}bpBin.${nPtsGenebody}PtsGenebody\", xlab=\"metagene\", ylab=\"Average pileups\" );
lines ( c((1-$nPts_uptss):(length(avg_bkt)-$nPts_uptss)), -avg_bkt_minus );
#lines( c($nPts_uptss,$nPts_uptss),range(avg_bkt), type=\"l\",,lty=2, col=\"red\" );
#lines( c($nPts_uptss+$nPtsGenebody,$nPts_uptss+$nPtsGenebody), range(avg_bkt), type=\"l\",lty=2, col=\"red\" );
dev.off()
""" >tmp.aggplot.metagenebucket.R

Rscript tmp.aggplot.metagenebucket.R ${dir}/${inprefix}.plus.pirna.genebucket ${dir}/${inprefix}.minus.pirna.genebucket ${dir}/output/$inprefix

#rm ${inprefix}.bed.pirna
#rm ${inprefix}.bed.pirna.clip.bed
#rm ${inprefix}.pirna.norm.bedGraph
