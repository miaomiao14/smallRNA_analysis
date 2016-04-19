#!/bin/bash
# pipeline for smallRNAs
# Author Jia
# March 2009
#modified by Wei 2013

#1.extract the inserts
insertsname=${1##*/}
dir=${1%/$insertsname}

inserts=$1
adaptor=$2
Min=$3
Max=$4
mm=$5
F=$6
EMAIL=$7
ID=$8

cp $inserts $inserts.uniq.reads

#1. map to the virus
echo "mapping to the virus..." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
run_bowtie.pl $inserts.uniq.reads $mm dmel_virus virus
virus_table.pl $inserts.virus > $dir/output/$insertsname.virus_table
exmatch.pl $inserts.virus $inserts.uniq.reads > $inserts.temp
mv $inserts.temp $inserts.uniq.reads
echo `date` >> $dir/output/$insertsname.log
echo "mapping to the virus done" >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log


#2.map to the genome 
echo "mapping to the genome...."  >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
if [ "$mm" == '0' ];
then
run_bowtie.pl $inserts.uniq.reads 1 dmel_TAS match2_all.out
run_bowtie_mismatchstats.pl $inserts.uniq.reads > $dir/output/$insertsname.mismatchstats.table
run_bowtie.pl $inserts.uniq.reads 0 dmel_TAS match2_all.out
else
run_bowtie.pl $inserts.uniq.reads 1 dmel_TAS match2_all.out
run_bowtie_mismatchstats.pl $inserts.uniq.reads > $dir/output/$insertsname.mismatchstats.table
fi
uniqmap.pl $inserts.match2_all.out > $inserts.uniqmap.match2_all.out
echo `date` >> $dir/output/$insertsname.log
echo "mapping to the genome done" >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log

if [ "$F" == 'Y' ];
then 
 grep -v YHet $inserts.match2_all.out > $inserts.xY.match2_all.out
 mv $inserts.xY.match2_all.out $inserts.match2_all.out
 grep -v YHet $inserts.uniqmap.match2_all.out > $inserts.uniqmap.xY.match2_all.out
 mv $inserts.uniqmap.xY.match2_all.out $inserts.uniqmap.match2_all.out
fi
cut -f1,2 $inserts.match2_all.out  |uniq.lines+ 0 > $inserts.match2_all.out.uniq.reads
cut -f1,2 $inserts.uniqmap.match2_all.out  |uniq.lines+ 0 > $inserts.uniqmap.match2_all.out.uniq.reads

# make bam 
#insert2bam.sh $inserts.uniq.reads

#3. map to transposon
echo "map to the transposons...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
matchall2normbed $inserts.match2_all.out $inserts.norm.bed
normbed2mapper /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/transposons+mst40+suffix+dodeca+stellate+TARTA+P.map $inserts.norm.bed map> $inserts.transposon.mapper
mapper2mapper2+ $inserts.transposon.mapper >  $inserts.transposon.mapper.temp
mod_stellate.pl  $inserts.transposon.mapper.temp > $inserts.transposon.mapper2
rm $inserts.transposon.mapper.temp
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.transposon.mapper2 > $inserts.uniqmap.transposon.mapper2
transposonlist.pl $inserts.transposon.mapper2 > $dir/output/$insertsname.transposon.list
transposonlist.pl $inserts.uniqmap.transposon.mapper2 >  $dir/output/$insertsname.uniqmap.transposon.list
cut -f1,2 $inserts.transposon.mapper2 | uniq.lines+ 0 > $inserts.transposon.mapper2.uniq.reads
cut -f1,2 $inserts.uniqmap.transposon.mapper2 | uniq.lines+ 0 > $inserts.uniqmap.transposon.mapper2.uniq.reads
echo `date` >> $dir/output/$insertsname.log
echo "map to the transposons done"  >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log

#4. map to everything else
echo "map to ncRNAs,miRNA,intron,exons,intergenic...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
#mapping.sh $inserts.uniq.reads $mm
# knownRNAs
echo "map to knownRNAs...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
maptocat $inserts.uniq.reads /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/knownRNAs.cat /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/knownRNAs.len $inserts.knownRNA.mapper
match.pl $inserts.uniqmap.match2_all.out.uniq.reads  $inserts.knownRNA.mapper > $inserts.uniqmap.knownRNA.mapper
echo `date` >> $dir/output/$insertsname.log
echo "map to knownRNAs done" >> $dir/output/$insertsname.log
 
# hairpins
echo "map to hairpins...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
#TODO: update the hairpin annotation; use bowtie to do the mapping
maptocat $inserts.uniq.reads /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/hairpin.cat /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/hairpin.len $inserts.hairpin_all.mapper
grep dme $inserts.hairpin_all.mapper > $inserts.hairpin.mapper
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.hairpin.mapper > $inserts.uniqmap.hairpin.mapper
rm $inserts.hairpin_all.mapper
annotmiR.pl $inserts.hairpin.mapper > $inserts.annotmiR
echo `date` >> $dir/output/$insertsname.log
echo "map to hairpins done" >> $dir/output/$insertsname.log

# mRNAs
echo "map to mRNAs...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
#TODO: update the mRNA annotation; build bowtie indexes
run_bowtie.pl $inserts.uniq.reads $mm dmel_transcript mRNA
modmRNAannot.pl $inserts.mRNA > $inserts.mRNA.mapper
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.mRNA.mapper > $inserts.uniqmap.mRNA.mapper
mRNAlist.pl $inserts.mRNA.mapper > $dir/output/$insertsname.mRNA.list
mRNAlist.pl $inserts.uniqmap.mRNA.mapper > $dir/output/$insertsname.uniqmap.mRNA.list
rm $inserts.mRNA
echo `date` >> $dir/output/$insertsname.log
echo "map to mRNAs done" >> $dir/output/$insertsname.log

grep + $inserts.mRNA.mapper > $inserts.snmRNA.mapper
grep -v + $inserts.mRNA.mapper > $inserts.antimRNA.mapper
grep + $inserts.uniqmap.mRNA.mapper > $inserts.uniqmap.snmRNA.mapper
grep -v + $inserts.uniqmap.mRNA.mapper > $inserts.uniqmap.antimRNA.mapper

# exons
echo "map to exons...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
#TODO: update the exon annotation; build bowtie indexes
run_bowtie.pl $inserts.uniq.reads $mm dmel_exon exon
bowtieout2mapper.pl $inserts.exon /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/exon.map > $inserts.exon.mapper
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.exon.mapper > $inserts.uniqmap.exon.mapper
rm $inserts.exon
echo `date` >> $dir/output/$insertsname.log
echo "map to exons done" >> $dir/output/$insertsname.log

# introns
echo "map to introns...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
#TODO: update the intron annotation; build bowtie indexes
run_bowtie.pl $inserts.uniq.reads $mm dmel_intron intron
bowtieout2mapper.pl $inserts.intron /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/intron.map > $inserts.intron.mapper
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.intron.mapper > $inserts.uniqmap.intron.mapper
rm $inserts.intron
echo `date` >> $dir/output/$insertsname.log
echo "map to introns done" >> $dir/output/$insertsname.log

# intergenic
echo "map to intergenic...." >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
#TODO: update the intergenic annotation; build bowtie indexes
run_bowtie.pl $inserts.uniq.reads $mm dmel_intergenic intergenic
bowtieout2mapper.pl $inserts.intergenic /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/intergenic.map > $inserts.intergenic.mapper
match.pl $inserts.uniqmap.match2_all.out.uniq.reads $inserts.intergenic.mapper > $inserts.uniqmap.intergenic.mapper
rm $inserts.intergenic
echo `date` >> $dir/output/$insertsname.log
echo "map to intergenic done" >> $dir/output/$insertsname.log

rm $inserts.uniq.reads.bowtie.out
echo "map to ncRNAs,miRNA,intron,exons,intergenic done" 
echo >> $dir/output/$insertsname.log

echo "making stats and annotation tables...."  >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
paraFile=$dir/output/${RANDOM}.para
echo -ne "table.sh $inserts $dir/output &&" >>${paraFile}
echo -e "tablenorm.pl $inserts $dir/output" >>${paraFile}
echo -ne "table.sh $inserts.uniqmap $dir/output &&" >>${paraFile}
echo -e "tablenorm.pl $inserts.uniqmap $dir/output" >>${paraFile}
if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	#CPUN=`wc -l $paraFile |cut -f1 -d" "` && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands
fi
echo `date` >> $dir/output/$insertsname.log
echo "table finished" >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log

# endosiRNA
echo "endosiRNA" >> $dir/output/$insertsname.log
echo `date` >> $dir/output/$insertsname.log
paraFile=$dir/output/${RANDOM}.para
##11/14/2013
##these norm.bed are not true norm.bed files
echo -ne "normbed2mapper /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/cisNAT.map  $inserts.xkxh.norm.bed  map  >  $inserts.xkxh.norm.bed.cisNAT &&" >>${paraFile}
echo -e "lenselector $inserts.xkxh.norm.bed.cisNAT 21 >  $inserts.xkxh.norm.bed.cisNAT.seq.21nt" >>${paraFile}
echo -ne "normbed2mapper /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/structured_loci.map  $inserts.xkxh.norm.bed  map  > $inserts.xkxh.norm.bed.structured_loci &&" >>${paraFile} 
echo -e "lenselector $inserts.xkxh.norm.bed.transposons 21 > $inserts.xkxh.norm.bed.transposons.seq.21nt" >>${paraFile}
echo -ne "normbed2mapper /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/transposons.map $inserts.xkxh.norm.bed map > $inserts.xkxh.norm.bed.transposons &&" >>${paraFile}
echo -e "lenrangeselector $inserts.xkxh.norm.bed.structured_loci 21 23 >  $inserts.xkxh.norm.bed.structured_loci.seq.21-23nt" >>${paraFile}
if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	#CPUN=`wc -l $paraFile |cut -f1 -d" "` && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands
fi
echo `date` >> $dir/output/$insertsname.log
echo "table finished" >> $dir/output/$insertsname.log
echo >> $dir/output/$insertsname.log

#5. make lendis figures
echo "making figures" >> $dir/output/$insertsname.log
paraFile=$dir/output/${RANDOM}.para

echo -ne "lendis $inserts.match2_all.out.uniq.reads > $inserts.match2_all.out.uniq.lendis && " >>${paraFile}
echo -e "RRR /home/wengz/pipelines/smallRNApipeline/pipeline_dm/R.source plot_lendis $inserts.match2_all.out.uniq.lendis ">>${paraFile}

echo -ne "lendis $inserts.match2_all.out.uniq.reads r > $inserts.match2_all.out.reads.lendis && ">>${paraFile}
echo -e "RRR /home/wengz/pipelines/smallRNApipeline/pipeline_dm/R.source plot_lendis $inserts.match2_all.out.reads.lendis ">>${paraFile}

echo -ne "lendis $inserts.xk.match2_all.out.uniq.reads > $inserts.xk.match2_all.out.uniq.lendis && ">>${paraFile}
echo -e "RRR /home/wengz/pipelines/smallRNApipeline/pipeline_dm/R.source plot_lendis $inserts.xk.match2_all.out.uniq.lendis ">>${paraFile}

echo -ne "lendis $inserts.xk.match2_all.out.uniq.reads r > $inserts.xk.match2_all.out.reads.lendis && ">>${paraFile}
echo -e "RRR /home/wengz/pipelines/smallRNApipeline/pipeline_dm/R.source plot_lendis $inserts.xk.match2_all.out.reads.lendis ">>${paraFile}

echo -ne "lendis $inserts.xkxh.match2_all.out.uniq.reads > $inserts.xkxh.match2_all.out.uniq.lendis && ">>${paraFile}
echo -e "RRR /home/wengz/pipelines/smallRNApipeline/pipeline_dm/R.source plot_lendis $inserts.xkxh.match2_all.out.uniq.lendis ">>${paraFile}

echo -ne "lendis $inserts.xkxh.match2_all.out.uniq.reads r > $inserts.xkxh.match2_all.out.reads.lendis && ">>${paraFile}
echo -e "RRR /home/wengz/pipelines/smallRNApipeline/pipeline_dm/R.source plot_lendis $inserts.xkxh.match2_all.out.reads.lendis ">>${paraFile}

echo -e "lendis2.pl $inserts.xkxh.transposon.mapper2 > $inserts.xkxh.transposon.mapper2.lendis2 &&">>${paraFile}
echo -e "RRR /home/wengz/pipelines/smallRNApipeline/pipeline_dm/R.source plot_lendis2 $inserts.xkxh.transposon.mapper2.lendis2 ">>${paraFile}

if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	#CPUN=`wc -l $paraFile |cut -f1 -d" "` && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands
fi
#6. readhist
#readhist $inserts.xkxh.match2_all.uniq.reads > $inserts.xkxh.match2_all.uniq.reads.readhist
#RRR /home/xuj1/pipeline/R.source plot_readhist $inserts.xkxh.match2_all.uniq.reads.readhist

echo "DONE" >> $dir/output/$insertsname.log

mv $dir/*.pdf $dir/output/
mv $dir/*.lendis $dir/output/
mv $dir/*.lendis2 $dir/output/

mail -s "Your job_$ID has completed" $EMAIL < /home/wengz/pipelines/smallRNApipeline/pipeline_dm/common/job_complete_message
echo "mismatch table" > $dir/internal_message
cat $dir/output/$insertsname.mismatchstats.table >> $dir/internal_message
echo >> $dir/internal_message
echo "stats reads table" >> $dir/internal_message
cat $dir/output/$insertsname.annot_table.reads >> $dir/internal_message
echo >> $dir/internal_message
echo "stats species table" >> $dir/internal_message
cat $dir/output/$insertsname.annot_table.species >> $dir/internal_message
echo >> $dir/internal_message
echo "stats reads table normalized by NTM and NTA" >> $dir/internal_message
cat $dir/output/$insertsname.annot_table.norm.reads >> $dir/internal_message
echo >> $dir/internal_message
echo "stats species table normalized by NTM and NTA" >> $dir/internal_message
cat $dir/output/$insertsname.annot_table.norm.species >> $dir/internal_message
echo >> $dir/internal_message
mail -s "job_$ID finished" "zzpipeline.admin@gmail.com" < $dir/internal_message
rm $dir/internal_message
rm $inserts.exon.mapper 
rm $inserts.intron.mapper 
rm $inserts.intergenic.mapper 
rm $inserts.mRNA.mapper
rm $inserts.snmRNA.mapper
rm $inserts.antimRNA.mapper

rm $inserts.uniqmap.exon.mapper
rm $inserts.uniqmap.intron.mapper
rm $inserts.uniqmap.intergenic.mapper
rm $inserts.uniqmap.mRNA.mapper
rm $inserts.uniqmap.snmRNA.mapper
rm $inserts.uniqmap.antimRNA.mapper

echo `date` >> $dir/output/$insertsname.log
echo "gzip the intermediate results...."  >> $dir/output/$insertsname.log
paraFile=$dir/output/${RANDOM}.para
for i in `ls -l $dir | egrep -v '^d' | awk '{print $9}'`
do
	echo -e "gzip ${dir}/$i" >>${paraFile}
done
if [[ ! -f ${paraFile}.completed ]] || [[ -f $$paraFile.failed_commands ]]
then
	#CPUN=`wc -l $paraFile |cut -f1 -d" "` && \
	ParaFly -c $paraFile -CPU 8 -failed_cmds $paraFile.failed_commands
fi
echo `date` >> $dir/output/$insertsname.log
echo "gzip the intermediate results done...."  >> $dir/output/$insertsname.log
