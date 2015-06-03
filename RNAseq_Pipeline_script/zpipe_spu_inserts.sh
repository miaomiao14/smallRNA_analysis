#!/bin/bash
# pipeline for smallRNAs spu

#1.extract the inserts
insertsname=${1##*/}
dir=${1%/$insertsname}

inserts=$1
adaptor=$2
Min=$3
Max=$4
mm=$5
EMAIL=$6
ID=$7
tissue=$8

mv $inserts $inserts.uniq.reads

#2.map to the genome 
echo "mapping to the genome...."  >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
/home/wangw1/pipeline_spu/run_bowtie_spu_un.pl $inserts.uniq.reads $mm /home/wangw1/pipeline_spu/common/Spur_3.1.LinearScaffold match2_all.out
uniqmap.pl $inserts.match2_all.out > $inserts.uniqmap.match2_all.out
echo `date` >> $dir/output/$insertsname.log
echo "mapping to the genome done" >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log

cut -f1,2 $inserts.match2_all.out  |uniq.lines+ 0 > $inserts.match2_all.out.uniq.reads
cut -f1,2 $inserts.uniqmap.match2_all.out  |uniq.lines+ 0 > $inserts.uniqmap.match2_all.out.uniq.reads

##map to the exon-exon splicing junction use all raw reads
echo "raw reads mapping to the splicing junction..." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
/home/wangw1/pipeline_spu/run_bowtie_spu.pl $inserts.uniq.reads $mm /home/wangw1/pipeline_spu/common/spu_${tissue}_splicing_junction match2_sj.out
FLANKING=30
awk -v F=$FLANKING '{OFS="\t"; if($5>= (F-length($1)+2) && $5 <=30 ) print $0;}' $inserts.match2_sj.out > $inserts.match2_flank30sj.out && \
rm $inserts.match2_sj.out
echo "raw reads mapping to the splicing junction done..." >> $dir/output/$insertsname.log
cut -f1,2 $inserts.match2_flank30sj.out |uniq.lines+ 0 > $inserts.match2_flank30sj.out.uniq.reads
uniqmap.pl $inserts.match2_flank30sj.out > $inserts.uniqmap.match2_flank30sj.out
cut -f1,2 $inserts.uniqmap.match2_flank30sj.out |uniq.lines+ 0 > $inserts.uniqmap.match2_flank30sj.out.uniq.reads

##TODO: map to the exon-intron splicing junctions use all raw reads

##TODO: map to the intron-exon splicing junctions use all raw reads


#3. map to transposon
echo "map to the transposons...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
matchall2normbed $inserts.match2_all.out $inserts.norm.bed
sed 's/Scaffold/chr/g' $inserts.norm.bed > $inserts.norm.bed.temp
normbed2mapper /home/wangw1/pipeline_spu/common/Baylor-Repeats.spurep_18.02.uniq.sort.grouped.map $inserts.norm.bed.temp map | sed 's/chr/Scaffold/g' > $inserts.transposon.mapper
mapper2mapper2+ $inserts.transposon.mapper > $inserts.transposon.mapper2
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.transposon.mapper2 > $inserts.uniqmap.transposon.mapper2
transposonlist.pl $inserts.transposon.mapper2 > $dir/output/$insertsname.transposon.list
transposonlist.pl $inserts.uniqmap.transposon.mapper2 > $dir/output/$insertsname.uniqmap.transposon.list
cut -f1,2 $inserts.transposon.mapper2 | uniq.lines+ 0 > $inserts.transposon.mapper2.uniq.reads
cut -f1,2 $inserts.uniqmap.transposon.mapper2 | uniq.lines+ 0 > $inserts.uniqmap.transposon.mapper2.uniq.reads
echo `date` >> $dir/output/$insertsname.log
echo "map to the transposons done"  >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log

#4. map to everything else
echo "map to ncRNAs,miRNA,Qiangintron,Qiangexons,intergenic...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
#mapping.sh $inserts.uniq.reads $mm
# knownRNAs
echo "map to knownRNAs...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
match.pl /home/xuj1/pipeline_spu/common/ncRNA_mapping.seq $inserts.uniq.reads > $inserts.knownRNA.mapper
match.pl /home/xuj1/pipeline_spu/common/ncRNA.mapper $inserts.uniqmap.match2_all.out.uniq.reads > $inserts.uniqmap.knownRNA.mapper
echo `date` >> $dir/output/$insertsname.log
echo "map to knownRNAs done" >> $dir/output/$insertsname.log

# hairpins
echo "map to hairpins...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log

/home/wangw1/pipeline_spu/run_bowtie_spu.pl $inserts.uniq.reads $mm /home/wangw1/indexes/Spu_hairpin hairpin.mapper 
grep "+" $inserts.hairpin.mapper > $inserts.hairpin.mapper.temp
mv $inserts.hairpin.mapper.temp $inserts.hairpin.mapper
cut -f1,2 $inserts.hairpin.mapper | uniq.lines+ 0 > $inserts.hairpin.mapper.uniq.reads

awk '{print $1,$2,$5,$6,$3}' $inserts.hairpin.mapper | tr -s ' ' '\t' > $inserts.hairpin.mapper.mod;
/home/wangw1/pipeline_spu/annotmiR_spu.pl /home/xuj1/nearline/spu_miRNA/all_miRNA.fa $inserts.hairpin.mapper.mod | uniq.lines+ 0 > $inserts.hairpin.mapper.annotmiR
miRlist.pl $inserts.hairpin.mapper.annotmiR > $dir/output/$insertsname.hairpin.mapper.annotmiR.list
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.hairpin.mapper > $inserts.uniqmap.hairpin.mapper

#map to different classes
/home/xuj1/pipeline_spu/run_bowtie.pl $inserts.uniq.reads $mm /home/wangw1/indexes/Spu_hairpin.c1 c1.hairpin.mapper
/home/xuj1/pipeline_spu/run_bowtie.pl $inserts.uniq.reads $mm /home/wangw1/indexes/Spu_hairpin.c2 c2.hairpin.mapper
/home/xuj1/pipeline_spu/run_bowtie.pl $inserts.uniq.reads $mm /home/wangw1/indexes/Spu_hairpin.c3 c3.hairpin.mapper
/home/xuj1/pipeline_spu/run_bowtie.pl $inserts.uniq.reads $mm /home/wangw1/indexes/Spu_hairpin.c4 c4.hairpin.mapper
grep "+" $inserts.c1.hairpin.mapper | awk '{print $1,$2,$5,$6,$3}' | tr -s ' ' '\t' > $inserts.c1.hairpin.mapper.mod
grep "+" $inserts.c2.hairpin.mapper | awk '{print $1,$2,$5,$6,$3}' | tr -s ' ' '\t' > $inserts.c2.hairpin.mapper.mod
grep "+" $inserts.c3.hairpin.mapper | awk '{print $1,$2,$5,$6,$3}' | tr -s ' ' '\t' > $inserts.c3.hairpin.mapper.mod
grep "+" $inserts.c4.hairpin.mapper | awk '{print $1,$2,$5,$6,$3}' | tr -s ' ' '\t' > $inserts.c4.hairpin.mapper.mod
/home/wangw1/pipeline_spu/annotmiR_spu.pl /home/xuj1/nearline/spu_miRNA/c1.miRNA.fa $inserts.c1.hairpin.mapper.mod | uniq.lines+ 0 > $inserts.c1.annotmiR
/home/wangw1/pipeline_spu/annotmiR_spu.pl /home/xuj1/nearline/spu_miRNA/c2.miRNA.fa $inserts.c2.hairpin.mapper.mod | uniq.lines+ 0 > $inserts.c2.annotmiR
/home/wangw1/pipeline_spu/annotmiR_spu.pl /home/xuj1/nearline/spu_miRNA/c3.miRNA.fa $inserts.c3.hairpin.mapper.mod | uniq.lines+ 0 > $inserts.c3.annotmiR
/home/wangw1/pipeline_spu/annotmiR_spu.pl /home/xuj1/nearline/spu_miRNA/c4.miRNA.fa $inserts.c4.hairpin.mapper.mod | uniq.lines+ 0 > $inserts.c4.annotmiR

echo `date` >> $dir/output/$insertsname.log
echo "map to hairpins done" >> $dir/output/$insertsname.log

# mRNAs
echo "map to mRNAs...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
#/home/xuj1/pipeline_spu/run_bowtie.pl $inserts.uniq.reads $mm Spu_mRNA_NCBI mRNA.mapper
#match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.mRNA.mapper > $inserts.uniqmap.mRNA.mapper

####
##WEI modify here 02/27
awk 'BEGIN{OFS="\t"}/track/{next};{print $1,$2-1,$3,$6,$7,$4,$5}' $inserts.norm.bed >$inserts.bed
bedtools intersect -a $inserts.bed -b /home/wangw1/pipeline_spu/common/Transcriptome_revised.gene.gff3 -f 1.0 -wb > $inserts.Qianggene.mapper.temp
awk 'BEGIN{OFS="\t"} { if($6 == $14){d="sense";} else {d="antisense";} split($16,a,";");gsub(/ID="/,"",a[1]);gsub(/"/,"",a[1]); print $7,$4,$1":"$2+1"-"$3"("$6")",d,$16,a[1],$5}' $inserts.Qianggene.mapper.temp > $inserts.Qianggene.mapper
rm $inserts.Qianggene.mapper.temp

mapper2mapper2+ $inserts.Qianggene.mapper > $inserts.Qianggene.mapper2
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.Qianggene.mapper2 > $inserts.uniqmap.Qianggene.mapper2

genelist.pl  $inserts.Qianggene.mapper2 >  $dir/output/$insertsname.Qianggene.list
cut -f1,2 $inserts.Qianggene.mapper2 | uniq.lines+ 0 > $inserts.Qianggene.mapper2.uniq.reads

genelist.pl $inserts.uniqmap.Qianggene.mapper2 > $dir/output/$insertsname.uniqmap.Qianggene.list
cut -f1,2 $inserts.uniqmap.Qianggene.mapper2 | uniq.lines+ 0 > $inserts.uniqmap.Qianggene.mapper2.uniq.reads


#intersectBed can intersect gtf file format; which incorporate the expression level of transcripts
bedtools intersect -a $inserts.bed -b /home/wangw1/pipeline_spu/common/spu_${tissue}_after_star_xrRNA_withSJDB_tophat_com_STAR_withRescued_cufflinks_o100_QiangGTF_unique.abundantisoform.sorted.gtf -f 1.0 -wb > $inserts.STARCUFFgene.mapper.temp
awk 'BEGIN{FS="\t";OFS="\t"} { if($6 == $14){d="sense";} else {d="antisense";} split($16,a,";");gsub(/transcript_id "/,"",a[2]);gsub(/"/,"",a[2]); print $7,$4,$1":"$2+1"-"$3"("$6")",d,$16,a[2],$5}' $inserts.STARCUFFgene.mapper.temp > $inserts.STARCUFFgene.mapper
rm $inserts.STARCUFFgene.mapper.temp
mapper2mapper2+ $inserts.STARCUFFgene.mapper >$inserts.STARCUFFgene.mapper2
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.STARCUFFgene.mapper2 >$inserts.uniqmap.STARCUFFgene.mapper2
genelistnew.pl $inserts.STARCUFFgene.mapper2 >$dir/output/$insertsname.STARCUFFgene.list
cut -f1,2 $inserts.STARCUFFgene.mapper2 | uniq.lines+ 0 >$inserts.STARCUFFgene.mapper2.uniq.reads
genelistnew.pl $inserts.uniqmap.STARCUFFgene.mapper2 >$dir/output/$insertsname.uniqmap.STARCUFFgene.list
cut -f1,2 $inserts.uniqmap.STARCUFFgene.mapper2 | uniq.lines+ 0 > $inserts.uniqmap.STARCUFFgene.mapper2.uniq.reads


# strand specific intersect
bedtools intersect -a $inserts.bed -b /home/wangw1/pipeline_spu/common/spu_${tissue}_after_star_xrRNA_withSJDB_tophat_com_STAR_withRescued_cufflinks_o100_QiangGTF_unique.abundantisoform.sorted.gtf -f 1.0 -wb -s > $inserts.STARCUFFgene.mapper.sense.temp
awk 'BEGIN{FS="\t";OFS="\t"} {split($16,a,";");gsub(/transcript_id "/,"",a[2]);gsub(/"/,"",a[2]); print $7,$4,$1":"$2+1"-"$3"("$6")","sense",$16,a[2],$5}' $inserts.STARCUFFgene.mapper.sense.temp >$inserts.STARCUFFgene.mapper.sense
rm $inserts.STARCUFFgene.mapper.sense.temp
bedtools intersect -a $inserts.bed -b /home/wangw1/pipeline_spu/common/spu_${tissue}_after_star_xrRNA_withSJDB_tophat_com_STAR_withRescued_cufflinks_o100_QiangGTF_unique.abundantisoform.sorted.gtf -f 1.0 -wb -S > $inserts.STARCUFFgene.mapper.antisense.temp
awk 'BEGIN{FS="\t";OFS="\t"} {split($16,a,";");gsub(/transcript_id "/,"",a[2]);gsub(/"/,"",a[2]); print $7,$4,$1":"$2+1"-"$3"("$6")","antisense",$16,a[2],$5}' $inserts.STARCUFFgene.mapper.antisense.temp >$inserts.STARCUFFgene.mapper.antisense
rm $inserts.STARCUFFgene.mapper.antisense.temp
cat $inserts.STARCUFFgene.mapper.sense $inserts.STARCUFFgene.mapper.antisense >$inserts.STARCUFFgene.strandspecific.mapper
rm $inserts.STARCUFFgene.mapper.sense $inserts.STARCUFFgene.mapper.antisense
mapper2mapper2+ $inserts.STARCUFFgene.strandspecific.mapper >$inserts.STARCUFFgene.strandspecific.mapper2
match.pl $inserts.uniqmap.match2_all.out.uniq.reads  $inserts.STARCUFFgene.strandspecific.mapper2 >$inserts.uniqmap.STARCUFFgene.strandspecific.mapper2
/home/wangw1/pipeline_spu/genelistnew.pl $inserts.STARCUFFgene.strandspecific.mapper2 >$dir/output/$insertsname.STARCUFFgene.strandspecific.mapper2.genelist
/home/wangw1/pipeline_spu/genelistnew.pl $inserts.uniqmap.STARCUFFgene.strandspecific.mapper2 >$dir/output/$insertsname.uniqmap.STARCUFFgene.strandspecific.mapper2 

bedtools intersect -a $inserts.bed -b /home/wangw1/pipeline_spu/common/spu_${tissue}_after_star_xrRNA_withSJDB_tophat_com_STAR_withRescued_cufflinks_o100_QiangGTF_unique.abundantisoformwithexon.exon.gtf -f 1.0 -wb > $inserts.STARCUFFexon.mapper.temp
awk 'BEGIN{OFS="\t"} { if($6 == $14){d="sense";} else {d="antisense";} split($16,a,";");gsub(/"/,"",a[1]);gsub(/"/,"",a[1]); print $7,$4,$1":"$2+1"-"$3"("$6")",d,$16,a[1],$5}' $inserts.STARCUFFexon.mapper.temp > $inserts.STARCUFFexon.mapper
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.STARCUFFexon.mapper >$inserts.uniqmap.STARCUFFexon.mapper

bedtools intersect -a $inserts.bed -b /home/wangw1/pipeline_spu/common/spu_${tissue}_after_star_xrRNA_withSJDB_tophat_com_STAR_withRescued_cufflinks_o100_QiangGTF_unique.abundantisoformwithexon.intron.bed -f 1.0 -wb > $inserts.STARCUFFintron.mapper.temp
awk 'BEGIN{OFS="\t"} { if($6 == $12){d="sense";} else {d="antisense";} print $7,$4,$1":"$2+1"-"$3"("$6")",d,$10,$10,$5}' $inserts.STARCUFFintron.mapper.temp >$inserts.STARCUFFintron.mapper
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.STARCUFFintron.mapper >$inserts.uniqmap.STARCUFFintron.mapper


intersectBed -a $inserts.bed -b /home/wangw1/pipeline_spu/common/Transcriptome_revised.exon.sort.grouped.gff3 -f 1.0 -wb > $inserts.Qiangexon.mapper.temp
awk 'BEGIN{OFS="\t"} { if($6 == $14){d="sense";} else {d="antisense";} split($16,a,";");gsub(/ID="/,"",a[1]);gsub(/"/,"",a[1]); print $7,$4,$1":"$2+1"-"$3"("$6")",d,$16,a[1],$5}' $inserts.Qiangexon.mapper.temp > $inserts.Qiangexon.mapper
rm $inserts.Qiangexon.mapper.temp
#mapper2mapper2+ $inserts.Qiangexon.mapper > $inserts.Qiangexon.mapper2
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.Qiangexon.mapper > $inserts.uniqmap.Qiangexon.mapper



intersectBed -a $inserts.bed -b /home/wangw1/pipeline_spu/common/Transcriptome_revised.intron.sort.grouped.gff3 -f 1.0 -wb > $inserts.Qiangintron.mapper.temp
awk 'BEGIN{OFS="\t"} { if($6 == $14){d="sense";} else {d="antisense";} split($16,a,";");gsub(/ID="/,"",a[1]);gsub(/"/,"",a[1]);print $7,$4,$1":"$2+1"-"$3"("$6")",d,$16,a[1],$5}' $inserts.Qiangintron.mapper.temp > $inserts.Qiangintron.mapper
rm $inserts.Qiangintron.mapper.temp
#mapper2mapper2+ $inserts.Qiangintron.mapper >$inserts.Qiangintron.mapper2
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.Qiangintron.mapper >$inserts.uniqmap.Qiangintron.mapper

rm $inserts.bed

normbed2mapper /home/xuj1/pipeline_spu/common/Spu_GLEAN_3UTR.map $inserts.norm.bed.temp map | sed 's/chr/Scaffold/g' > $inserts.3UTR.mapper
rm $inserts.norm.bed.temp
#mapper2mapper2+ $inserts.3UTR.mapper >$inserts.3UTR.mapper2
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.3UTR.mapper > $inserts.uniqmap.3UTR.mapper
#mapper2mapper2+ $inserts.uniqmap.3UTR.mapper >$inserts.uniqmap.3UTR.mapper2

#map to splicing junctions#

echo `date` >> $dir/output/$insertsname.log
echo "map to mRNAs done" >> $dir/output/$insertsname.log






echo `date` >> $dir/output/$insertsname.log
echo "map non-genome mapping reads to new transcriptome assembled by STAR and cufflinks..." >> $dir/output/$insertsname.log
exmatch.pl $inserts.hairpin.mapper $inserts.um.reads >$inserts.um.xh.reads
match.pl $inserts.um.xh.reads $inserts.uniq.reads > $inserts.um.uniq.reads
rm $inserts.um.reads
/home/wangw1/pipeline_spu/run_bowtie_spu_un.pl $inserts.um.uniq.reads $mm /home/wangw1/pipeline_spu/common/spu_${tissue}_transcriptome match2_splicing.out
uniqmap.pl $inserts.um.match2_splicing.out >$inserts.uniqmap.um.match2_splicing.out
cat $inserts.um.match2_splicing.out |awk '{OFS="\t"}{print "transcript"$3,$5,$6,$1,$2,$4}' | bedtools sort -i stdin > $inserts.um.match2cufftranscriptome.out.bed
cat $inserts.uniqmap.um.match2_splicing.out |awk '{OFS="\t"}{print "transcript"$3,$5,$6,$1,$2,$4}' | bedtools sort -i stdin > $inserts.uniqmap.um.match2cufftranscriptome.out.bed

echo `date` >> $dir/output/$insertsname.log
match.pl $inserts.um.um.reads $inserts.uniq.reads > $inserts.um.um.uniq.reads
rm $inserts.um.um.reads
echo "map non-genome mapping reads to new transcriptome assembled by trinity..." >> $dir/output/$insertsname.log
/home/wangw1/pipeline_spu/run_bowtie_spu_un.pl $inserts.um.um.uniq.reads $mm /home/wangw1/pipeline_spu/common/spu_${tissue}_transcriptome_trinity match2_splicing.out 
uniqmap.pl $inserts.um.um.match2_splicing.out > $inserts.uniqmap.um.um.match2_splicing.out
cat $inserts.um.um.match2_splicing.out |awk '{OFS="\t"}{print $3,$5,$6,$1,$2,$4}' |bedtools sort -i stdin  >$inserts.um.match2trinitytranscriptome.out.bed
cat $inserts.uniqmap.um.um.match2_splicing.out |awk '{OFS="\t"}{print "transcript"$3,$5,$6,$1,$2,$4}' | bedtools sort -i stdin > $inserts.uniqmap.um.match2trinitytranscriptome.out.bed
#echo `date` >> $dir/output/$insertsname.log
#echo "concatate rescue reads and convert it to bed format.." >> $dir/output/$insertsname.log
#cat $inserts.um.match2cufftranscriptome.out.bed $inserts.um.match2trinitytranscriptome.out.bed > $inserts.um.match2transcriptome.out.bed
#cat $inserts.um.match2transcriptome.out.bed |sort -k1,1 -k2,2n -k3,3n -k6,6 | bedtools groupby -g 1,2,3,4,6 -c 5 -o count -full > $inserts.um.match2transcriptome.out.merged.bed
#awk '{OFS="\t"} { if($7==1) print $0}' $inserts.um.match2transcriptome.out.bed > $inserts.um.match2transcriptome.uniqmap.bed  

echo `date` >> $dir/output/$insertsname.log
echo "map non-genome mapping reads to new transcriptome assembled by trinity done..." >> $dir/output/$insertsname.log

echo "making stats and annotation tables...."  >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
/home/wangw1/pipeline_spu/table_spu.sh $inserts
/home/wangw1/pipeline_spu/table_spu.sh $inserts.uniqmap
mv $inserts*table $dir/output
mv $inserts*stats* $dir/output
echo `date` >> $dir/output/$insertsname.log
echo "table finished" >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log

#5. make lendis figures
echo "making figures" >> $dir/output/$insertsname.log
lendis $inserts.match2_all.out.uniq.reads > $inserts.match2_all.out.uniq.lendis
lendis $inserts.match2_all.out.uniq.reads r > $inserts.match2_all.out.reads.lendis
RRR /home/xuj1/pipeline/R.source plot_lendis $inserts.match2_all.out.uniq.lendis
RRR /home/xuj1/pipeline/R.source plot_lendis $inserts.match2_all.out.reads.lendis

lendis $inserts.xk.match2_all.out.uniq.reads > $inserts.xk.match2_all.out.uniq.lendis
lendis $inserts.xk.match2_all.out.uniq.reads r > $inserts.xk.match2_all.out.reads.lendis
RRR /home/xuj1/pipeline/R.source plot_lendis $inserts.xk.match2_all.out.uniq.lendis
RRR /home/xuj1/pipeline/R.source plot_lendis $inserts.xk.match2_all.out.reads.lendis

lendis $inserts.xkxh.match2_all.out.uniq.reads > $inserts.xkxh.match2_all.out.uniq.lendis
lendis $inserts.xkxh.match2_all.out.uniq.reads r > $inserts.xkxh.match2_all.out.reads.lendis
RRR /home/xuj1/pipeline/R.source plot_lendis $inserts.xkxh.match2_all.out.uniq.lendis
RRR /home/xuj1/pipeline/R.source plot_lendis $inserts.xkxh.match2_all.out.reads.lendis

/home/xuj1/pipeline_spu/lendis2.pl $inserts.xkxh.transposon.mapper2 > $inserts.xkxh.transposon.mapper2.lendis2
RRR /home/xuj1/pipeline/R.source plot_lendis2 $inserts.xkxh.transposon.mapper2.lendis2
/home/xuj1/pipeline_spu/lendis2.pl $inserts.xkxh.Qianggene.mapper2 > $inserts.xkxh.Qianggene.mapper2.lendis2
RRR /home/xuj1/pipeline/R.source plot_lendis2 $inserts.xkxh.Qianggene.mapper2.lendis2

/home/xuj1/pipeline_spu/lendis2.pl $inserts.xkxh.STARCUFFgene.mapper2 > $inserts.xkxh.STARCUFFgene.mapper2.lendis2
RRR /home/xuj1/pipeline/R.source plot_lendis2 $inserts.xkxh.STARCUFFgene.mapper2.lendis2


echo "DONE" >> $dir/output/$insertsname.log

mv $dir/*.pdf $dir/output/
mv $dir/*.lendis $dir/output/
mv $dir/*.lendis2 $dir/output/

#rm $inserts.Qianggene.mapper
#rm $inserts.STARCUFFgene.mapper
rm $inserts.Qiangintron.mapper
rm $inserts.Qiangexon.mapper 
rm $inserts.3UTR.mapper
#rm $inserts.mRNA.mapper

rm $inserts.uniqmap.Qianggene.mapper
rm $inserts.uniqmap.Qiangintron.mapper
rm $inserts.uniqmap.Qiangexon.mapper
rm $inserts.uniqmap.3UTR.mapper
#rm $inserts.uniqmap.mRNA.mapper
echo `date` >> $dir/output/$insertsname.log
echo "metagene analysis...."  >> $dir/output/$insertsname.log
/home/wangw1/bin/normbed2metageneAgg.sh $inserts.xkxh.norm.bed $tissue
echo `date` >> $dir/output/$insertsname.log
echo "metagene analysis done...."  >> $dir/output/$insertsname.log

echo `date` >> $dir/output/$insertsname.log
echo "gzip the intermediate results...."  >> $dir/output/$insertsname.log 
for i in `ls -l $dir | egrep -v '^d' | awk '{print $9}'`
do 
  gzip ${dir}/$i
done
echo `date` >> $dir/output/$insertsname.log
echo "gzip the intermediate results done...."  >> $dir/output/$insertsname.log