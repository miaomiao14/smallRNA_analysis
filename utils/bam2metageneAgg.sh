#!/bin/sh
# calculate the aggregation on "metagenes" (average over all genes)
# 2013-02-27, tushikui@gmail

inbam=$1      # Input bam file (single)
in=${1##*/}
dir=${1%/$in}
ouname=`basename $in .bam`
tissue=$2     # name for outputs

chrom_sizes=/home/wangw1/pipeline_spu/common/Spur_3.1.scaffoldsize # chrom sizes file (two col: chrom size)
gdb=spu
genebed=/home/wangw1/pipeline_spu/common/spu_${tissue}_after_star_xrRNA_withSJDB_tophat_com_STAR_withRescued_cufflinks_o100_QiangGTF_unique.sort.bed     # dmel-all-r5.45.gff.FlyBase.gene.bed
normFac=`grep ${inprefix} /home/wangw1/pipeline_spu/common/spu.norm.factor | cut -f2`     # *1000000/normFac (from samp.norm.factor.allmap)

tssUpLen=1000     # for "+": [TSS-1000,TSS]; "-":[TSS,TSS+1000]
ttsDnLen=1000     # for "+": [TTS,TTS+1000]; "-":[TTS-1000,TTS]
nPtsGenebody=1000   # num of points (intervals) in [TSS,TTS]
grid=20              # grid width for up/down-stream
d=150      # bp for shift or extension (for ChIP-seq reads)
#gdb="dm3"  # genome database

declare -i nPts_uptss=$tssUpLen/$grid
declare -i nPts_dntts=$ttsDnLen/$grid


## bam -> bed  (shift d only for ChIP-seq)
[ ! -f ${dir}/tmp.bed ] && \
bedtools bamtobed -i $inbam | awk -F "\t" '{ OFS="\t"; if($6=="-") { print $1,$2,$3,".",0,"-" } else { print $1,$2,$3,".",0,"+" } }' | bedClip stdin $chrom_sizes ${dir}/tmp.bed


## pileup: bed -> bedGraph
[ ! -f ${dir}/$ouname.bedGraph ] && \
cut -f1-3 ${dir}/tmp.bed | sort -k1,1 -k2,2n | bedItemOverlapCount $gdb stdin -chromSize=$chrom_sizes | awk -F "\t" -v normFac=$normFac '{ OFS="\t"; print $1,$2,$3,$4*1000000/normFac; }' >${dir}/$ouname.bedGraph

## bedGraph -> bw
[ ! -f ${dir}/$ouname.bw ] && \
bedGraphToBigWig ${dir}/$ouname.bedGraph $chrom_sizes ${dir}/$ouname.bw


## get gene info (chr,tss,tts,name): reduce duplicate gene names
[ ! -f ${dir}/gene_tss_tts.bed ] && cut -f1-6 $genebed | sort -k4,4 -u | sort -k1,1 -k2,2n >${dir}/gene_tss_tts.bed

## get gene bucket value
fetchBWscore_genebucket.pl ${dir}/$ouname.bw ${dir}/gene_tss_tts.bed $tssUpLen $nPtsGenebody $ttsDnLen $grid >${dir}/$ouname.genebucket


## plot the metagene bucket value
echo """#!/home/tus/bin/Rscript
args <- commandArgs( TRUE );

data_file <- args[1];
ouname    <- args[2];

# read data: nGen*nWid
bkt_val <- read.table( data_file, sep=\"\t\", header=FALSE );

# average signals across rows (genes)
avg_bkt <- apply( bkt_val, 2, mean );

# plot the avg signals
pdf( paste( ouname, \"genebucket.pdf\", sep=\".\" ) );
plot( c((1-$nPts_uptss):(length(avg_bkt)-$nPts_uptss)), avg_bkt, type=\"l\", main=\"gene.aggplot.up$tssUpLen.down$ttsDnLen.${grid}bpBin.${nPtsGenebody}PtsGenebody\", xlab=\"metagene\", ylab=\"Average pileups\" );
#lines( c($nPts_uptss,$nPts_uptss),range(avg_bkt), type=\"l\", col=\"red\" );
#lines( c($nPts_uptss+$nPtsGenebody,$nPts_uptss+$nPtsGenebody), range(avg_bkt), type=\"l\", col=\"red\" );
dev.off()
""" >${dir}/tmp.aggplot.metagenebucket.R
Rscript ${dir}/tmp.aggplot.metagenebucket.R ${dir}/${ouname}.genebucket ${dir}/${ouname}




#if [ -f tmp_bkt.val ]; then
#	rm tmp_bkt.val
#fi
#cat gene_tss_tts.bed | while read chr chrStart chrEnd genename score strand
#do
#	if [ $strand == "+" ]; then
#		declare -i upchrStart=$chrStart-$tssUpLen
#		declare -i dnchrEnd=$chrEnd+$ttsDnLen
#		# upstream: [TSS-1000,TSS]
#		bigWigSummary ${ouname}.bw $chr $upchrStart $chrStart $nPts_uptss >t_uptss_bkt.val
#		# gene body: [TSS,TTS]
#		bigWigSummary ${ouname}.bw $chr $chrStart $chrEnd $nPtsGenebody >t_genebody_bkt.val
#		# downstream: [TTS,TTS+1000]
#		bigWigSummary ${ouname}.bw $chr $chrEnd $dnchrEnd $nPts_dntts >t_dntts_bkt.val
#
#		# join the three parts: tss+genebody+tts
#		paste t_uptss_bkt.val t_genebody_bkt.val t_dntss_bkt.val >> tmp_bkt.val
#	fi
#
#	if [ $strand == "-" ]; then
#		declare -i upchrEnd=$chrEnd+$tssUpLen
#		declare -i dnchrStart=$chrStart-$ttsDnLen
#		# upstream: [TSS,TSS+1000]
#		bigWigSummary ${ouname}.bw $chr $chrEnd $upchrEnd $nPts_uptss >t_uptss_bkt.val
#		# gene body: [TTS,TSS]
#		bigWigSummary ${ouname}.bw $chr $chrStart $chrEnd $nPtsGenebody >t_genebody_bkt.val
#		# downstream: [TTS-1000,TTS]
#		bigWigSummary ${ouname}.bw $chr $dnchrStart $chrStart $nPts_dntts >t_dntts_bkt.val
#
#		# join the three parts: tts+genebody+tss
#		paste t_dntss_bkt.val t_genebody_bkt.val t_uptss_bkt.val >tmp
#		
#		# reverse the value vector for '-' strand
#		awk -F "\t" '
#		{   for( i=NF; i>1; i-- ) 
#			{
#				printf( "%s\t", $i );
#			}
#			printf( "%s\n", $1 );
#		}' tmp >>tmp_bkt.val
#	fi
#done
#
## get the bucket value file
#mv tmp_bkt.val ${ouname}.genebucket
#cat tmp_bkt.val
