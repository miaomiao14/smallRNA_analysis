#!/bin/bash

#table.sh <insert>

insert_species=`wc -l $1.uniq.reads | cut -d " " -f1`
insert_reads=`sumcol $1.uniq.reads`

genome_species=`wc -l $1.match2_all.out.uniq.reads | cut -d " " -f1`
genome_reads=`sumcol $1.match2_all.out.uniq.reads`

exmatch.pl $1.knownRNA.mapper $1.match2_all.out.uniq.reads > $1.xk.match2_all.out.uniq.reads
exknown_species=`wc -l $1.xk.match2_all.out.uniq.reads | cut -d " " -f1`
exknown_reads=`sumcol $1.xk.match2_all.out.uniq.reads`


exmatch.pl $1.hairpin.mapper $1.xk.match2_all.out.uniq.reads > $1.xkxh.match2_all.out.uniq.reads
exkexh_species=`wc -l $1.xkxh.match2_all.out.uniq.reads | cut -d " " -f1`
exkexh_reads=`sumcol $1.xkxh.match2_all.out.uniq.reads`

let "nc_species=genome_species-exknown_species"
let "nc_reads=genome_reads-exknown_reads"

let "mir_species=exknown_species-exkexh_species"
let "mir_reads=exknown_reads-exkexh_reads"

match.pl $1.xkxh.match2_all.out.uniq.reads $1.match2_all.out > $1.xkxh.match2_all.out
matchall2normbed $1.xkxh.match2_all.out $1.xkxh.norm.bed
cut -f1,2 $1.xkxh.match2_all.out | uniq.lines+ 0 > $1.xkxh.match2_all.out.uniq.reads
match.pl $1.xkxh.match2_all.out.uniq.reads $1.transposon.mapper2 > $1.xkxh.transposon.mapper2
cut -f1,2 $1.xkxh.transposon.mapper2 | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.uniq.reads


#--------------PiRNA----------------------
lenrangeselector $1.xkxh.match2_all.out.uniq.reads 23 29 >  $1.xkxh.match2_all.out.23-29nt.uniq.reads
lenrangeselector $1.xkxh.transposon.mapper2 23 29 > $1.xkxh.transposon.mapper2.23-29nt
lenrangeselector $1.xkxh.transposon.mapper2.uniq.reads 23 29 > $1.xkxh.transposon.mapper2.23-29nt.uniq.reads

piRNA_species=`wc -l $1.xkxh.match2_all.out.23-29nt.uniq.reads| cut -d " " -f1`
piRNA_reads=`sumcol $1.xkxh.match2_all.out.23-29nt.uniq.reads`

grep antisense $1.xkxh.transposon.mapper2.23-29nt  | cut -f1,2 | uniq.lines+ 0 >  $1.xkxh.transposon.mapper2.23-29nt.antisense.uniq.reads
grep -v antisense $1.xkxh.transposon.mapper2.23-29nt | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.23-29nt.sense.uniq.reads

trans_piRNA_total_species=`wc -l $1.xkxh.transposon.mapper2.23-29nt.uniq.reads  | cut -d " " -f1`
trans_piRNA_total_reads=`sumcol $1.xkxh.transposon.mapper2.23-29nt.uniq.reads `

trans_piRNA_sense_species=`wc -l $1.xkxh.transposon.mapper2.23-29nt.sense.uniq.reads  | cut -d " " -f1`
trans_piRNA_sense_reads=`sumcol $1.xkxh.transposon.mapper2.23-29nt.sense.uniq.reads`

trans_piRNA_antisense_species=`wc -l $1.xkxh.transposon.mapper2.23-29nt.antisense.uniq.reads  | cut -d " " -f1`
trans_piRNA_antisense_reads=`sumcol $1.xkxh.transposon.mapper2.23-29nt.antisense.uniq.reads`

#------------------siRNA-------------------
lenselector $1.xkxh.match2_all.out.uniq.reads 21 > $1.xkxh.match2_all.out.21nt.uniq.reads
lenselector $1.xkxh.transposon.mapper2 21 > $1.xkxh.transposon.mapper2.21nt
lenselector $1.xkxh.transposon.mapper2.uniq.reads 21 > $1.xkxh.transposon.mapper2.21nt.uniq.reads

siRNA_species=`wc -l $1.xkxh.match2_all.out.21nt.uniq.reads| cut -d " " -f1`
siRNA_reads=`sumcol $1.xkxh.match2_all.out.21nt.uniq.reads`

grep antisense  $1.xkxh.transposon.mapper2.21nt | cut -f1,2 | uniq.lines+ 0 >  $1.xkxh.transposon.mapper2.21nt.antisense.uniq.reads
grep -v antisense  $1.xkxh.transposon.mapper2.21nt | cut -f1,2 | uniq.lines+ 0 >  $1.xkxh.transposon.mapper2.21nt.sense.uniq.reads

trans_siRNA_total_species=`wc -l $1.xkxh.transposon.mapper2.21nt.uniq.reads  | cut -d " " -f1`
trans_siRNA_total_reads=`sumcol $1.xkxh.transposon.mapper2.21nt.uniq.reads `

trans_siRNA_sense_species=`wc -l $1.xkxh.transposon.mapper2.21nt.sense.uniq.reads  | cut -d " " -f1`
trans_siRNA_sense_reads=`sumcol $1.xkxh.transposon.mapper2.21nt.sense.uniq.reads`

trans_siRNA_antisense_species=`wc -l $1.xkxh.transposon.mapper2.21nt.antisense.uniq.reads  | cut -d " " -f1`
trans_siRNA_antisense_reads=`sumcol $1.xkxh.transposon.mapper2.21nt.antisense.uniq.reads`

echo -e "inserts\tgenome_mapping\tncRNAs\texcluding_ncRNAs\tmiRNAs\texcluding_ncRNAs_and_hairpins\tpiRNAs\ttransposon_map_piRNAs\tsense\tantisense\tsiRNAs\ttransposon_map_siRNAs\tsense\tantisense" > $1_stats_table_species
echo -e "$insert_species\t$genome_species\t$nc_species\t$exknown_species\t$mir_species\t$exkexh_species\t$piRNA_species\t$trans_piRNA_total_species\t$trans_piRNA_sense_species\t$trans_piRNA_antisense_species\t$siRNA_species\t$trans_siRNA_total_species\t$trans_siRNA_sense_species\t$trans_siRNA_antisense_species" >> $1_stats_table_species
echo -e "inserts\tgenome_mapping\tncRNAs\texcluding_ncRNAs\tmiRNAs\texcluding_ncRNAs_and_hairpins\tpiRNAs\ttransposon_map_piRNAs\tsense\tantisense\tsiRNAs\ttransposon_map_siRNAs\tsense\tantisense" > $1_stats_table_reads
echo -e "$insert_reads\t$genome_reads\t$nc_reads\t$exknown_reads\t$mir_reads\t$exkexh_reads\t$piRNA_reads\t$trans_piRNA_total_reads\t$trans_piRNA_sense_reads\t$trans_piRNA_antisense_reads\t$siRNA_reads\t$trans_siRNA_total_reads\t$trans_siRNA_sense_reads\t$trans_siRNA_antisense_reads" >> $1_stats_table_reads

#---------------annotation---------------------------
exmatch2.pl $1.knownRNA.mapper $1.hairpin.mapper $1.snmRNA.mapper > $1.xkxh.snmRNA.mapper
exmatch2.pl $1.knownRNA.mapper $1.hairpin.mapper $1.antimRNA.mapper > $1.xkxh.antimRNA.mapper
exmatch2.pl $1.knownRNA.mapper $1.hairpin.mapper $1.exon.mapper > $1.xkxh.exon.mapper
exmatch2.pl $1.knownRNA.mapper $1.hairpin.mapper $1.intergenic.mapper > $1.xkxh.intergenic.mapper
exmatch2.pl $1.knownRNA.mapper $1.hairpin.mapper $1.intron.mapper > $1.xkxh.intron.mapper
exmatch2.pl $1.knownRNA.mapper $1.hairpin.mapper $1.mRNA.mapper > $1.xkxh.mRNA.mapper

cut -f1,2  $1.xkxh.snmRNA.mapper | uniq.lines+ 0 > $1.xkxh.snmRNA.mapper.uniq.reads
cut -f1,2  $1.xkxh.antimRNA.mapper | uniq.lines+ 0  > $1.xkxh.antimRNA.mapper.uniq.reads

match.pl $1.xkxh.transposon.mapper2.uniq.reads $1.xkxh.snmRNA.mapper.uniq.reads > $1.xkxh.snmRNA.overlaptrans.uniq.reads
match.pl $1.xkxh.transposon.mapper2.uniq.reads $1.xkxh.antimRNA.mapper.uniq.reads > $1.xkxh.antimRNA.overlaptrans.uniq.reads
match.pl $1.xkxh.match2_all.out.uniq.reads  $1.intergenic.mapper | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.match2_all.out.intergenic.uniq.reads
match.pl $1.xkxh.transposon.mapper2.uniq.reads $1.xkxh.match2_all.out.intergenic.uniq.reads > $1.xkxh.match2_all.out.intergenic.overlaptrans.uniq.reads
match.pl $1.xkxh.match2_all.out.uniq.reads  $1.intron.mapper | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.match2_all.out.intron.uniq.reads
match.pl $1.xkxh.transposon.mapper2.uniq.reads $1.xkxh.match2_all.out.intron.uniq.reads > $1.xkxh.match2_all.out.intron.overlaptrans.uniq.reads
match.pl $1.xkxh.match2_all.out.uniq.reads  $1.exon.mapper | cut -f1,2 | uniq.lines+ 0 >$1.xkxh.match2_all.out.exon.uniq.reads
match.pl $1.xkxh.transposon.mapper2.uniq.reads $1.xkxh.match2_all.out.exon.uniq.reads > $1.xkxh.match2_all.out.exon.overlaptrans.uniq.reads
exmatch5.pl $1.xkxh.snmRNA.mapper.uniq.reads $1.xkxh.antimRNA.mapper.uniq.reads $1.xkxh.match2_all.out.intergenic.uniq.reads $1.xkxh.match2_all.out.intron.uniq.reads $1.xkxh.match2_all.out.exon.uniq.reads $1.xkxh.match2_all.out.uniq.reads > $1.xkxh.match2_all.out.noannot.uniq.reads

grep Retroviral_elements $1.xkxh.transposon.mapper2 | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.LTR.uniq.reads
grep Foldback_elements $1.xkxh.transposon.mapper2 | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.Foldback.uniq.reads
grep IR-elements $1.xkxh.transposon.mapper2 | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.IR.uniq.reads
grep SINE-like_elements $1.xkxh.transposon.mapper2 | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.SINE.uniq.reads
grep non-LTR_retrotransposons  $1.xkxh.transposon.mapper2 | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.nonLTR.uniq.reads
grep repeats $1.xkxh.transposon.mapper2 | cut -f1,2 | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.repeats.uniq.reads

cat $1.xkxh.transposon.mapper2.SINE.uniq.reads $1.xkxh.transposon.mapper2.nonLTR.uniq.reads | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.non_LTR.uniq.reads
cat $1.xkxh.transposon.mapper2.Foldback.uniq.reads $1.xkxh.transposon.mapper2.IR.uniq.reads | uniq.lines+ 0 > $1.xkxh.transposon.mapper2.DNA.uniq.reads

LTR_species=`wc -l $1.xkxh.transposon.mapper2.LTR.uniq.reads | cut -d " " -f1`
LTR_reads=`sumcol $1.xkxh.transposon.mapper2.LTR.uniq.reads`

non_LTR_species=`wc -l $1.xkxh.transposon.mapper2.non_LTR.uniq.reads | cut -d " " -f1`
non_LTR_reads=`sumcol $1.xkxh.transposon.mapper2.non_LTR.uniq.reads`

DNA_species=`wc -l $1.xkxh.transposon.mapper2.DNA.uniq.reads | cut -d " " -f1`
DNA_reads=`sumcol $1.xkxh.transposon.mapper2.DNA.uniq.reads`

repeat_species=`wc -l $1.xkxh.transposon.mapper2.repeats.uniq.reads | cut -d " " -f1`
repeat_reads=`sumcol $1.xkxh.transposon.mapper2.repeats.uniq.reads`

sense_mRNA_species=`wc -l $1.xkxh.snmRNA.mapper.uniq.reads | cut -d " " -f1`
sense_mRNA_reads=`sumcol $1.xkxh.snmRNA.mapper.uniq.reads`

sense_mRNA_overlap_species=`wc -l $1.xkxh.snmRNA.overlaptrans.uniq.reads | cut -d " " -f1`
sense_mRNA_overlap_reads=`sumcol $1.xkxh.snmRNA.overlaptrans.uniq.reads`

let "sense_mRNA_nonoverlap_species=sense_mRNA_species-sense_mRNA_overlap_species"
let "sense_mRNA_nonoverlap_reads=sense_mRNA_reads-sense_mRNA_overlap_reads"

antisense_mRNA_species=`wc -l $1.xkxh.antimRNA.mapper.uniq.reads | cut -d " " -f1`
antisense_mRNA_reads=`sumcol $1.xkxh.antimRNA.mapper.uniq.reads`

antisense_mRNA_overlap_species=`wc -l $1.xkxh.antimRNA.overlaptrans.uniq.reads | cut -d " " -f1`
antisense_mRNA_overlap_reads=`sumcol $1.xkxh.antimRNA.overlaptrans.uniq.reads`

let "antisense_mRNA_nonoverlap_species=antisense_mRNA_species-antisense_mRNA_overlap_species"
let "antisense_mRNA_nonoverlap_reads=antisense_mRNA_reads-antisense_mRNA_overlap_reads"

intron_species=`wc -l $1.xkxh.match2_all.out.intron.uniq.reads | cut -d " " -f1`
intron_reads=`sumcol $1.xkxh.match2_all.out.intron.uniq.reads`

intron_overlap_species=`wc -l $1.xkxh.match2_all.out.intron.overlaptrans.uniq.reads | cut -d " " -f1`
intron_overlap_reads=`sumcol $1.xkxh.match2_all.out.intron.overlaptrans.uniq.reads`

let "intron_nonoverlap_species=intron_species-intron_overlap_species"
let "intron_nonoverlap_reads=intron_reads-intron_overlap_reads"

exon_species=`wc -l $1.xkxh.match2_all.out.exon.uniq.reads | cut -d " " -f1`
exon_reads=`sumcol $1.xkxh.match2_all.out.exon.uniq.reads`

exon_overlap_species=`wc -l $1.xkxh.match2_all.out.exon.overlaptrans.uniq.reads | cut -d " " -f1`
exon_overlap_reads=`sumcol $1.xkxh.match2_all.out.exon.overlaptrans.uniq.reads`

let "exon_nonoverlap_species=exon_species-exon_overlap_species"
let "exon_nonoverlap_reads=exon_reads-exon_overlap_reads"

intergenic_species=`wc -l $1.xkxh.match2_all.out.intergenic.uniq.reads | cut -d " " -f1`
intergenic_reads=`sumcol $1.xkxh.match2_all.out.intergenic.uniq.reads`

intergenic_overlap_species=`wc -l $1.xkxh.match2_all.out.intergenic.overlaptrans.uniq.reads | cut -d " " -f1`
intergenic_overlap_reads=`sumcol $1.xkxh.match2_all.out.intergenic.overlaptrans.uniq.reads`

let "intergenic_nonoverlap_species=intergenic_species-intergenic_overlap_species"
let "intergenic_nonoverlap_reads=intergenic_reads-intergenic_overlap_reads"


noannot_species=`wc -l $1.xkxh.match2_all.out.noannot.uniq.reads | cut -d " " -f1`
noannot_reads=`sumcol $1.xkxh.match2_all.out.noannot.uniq.reads`

echo -e "Annotation\tLTR\tnon_LTR\tDNA\trepeats\tsensemRNA\tsensemRNA_overlap_transposons\tsensemRNA_nonoverlap_transposons\tantisensemRNA\tantisensemRNA_overlap_transposons\tantisensemRNA_nonoverlap_transposons\tintron\tintron_overlap_transposons\tintron_nonoverlap_transposons\texon\texon_overlap_transposons\texon_nonoverlap_transposons\tintergenic\tintergenic_overlap_transposons\tintergenic_nonoverlap_transposons\tno annotation" > $1.annot_table.species
echo -e "All\t$LTR_species\t$non_LTR_species\t$DNA_species\t$repeat_species\t$sense_mRNA_species\t$sense_mRNA_overlap_species\t$sense_mRNA_nonoverlap_species\t$antisense_mRNA_species\t$antisense_mRNA_overlap_species\t$antisense_mRNA_nonoverlap_species\t$intron_species\t$intron_overlap_species\t$intron_nonoverlap_species\t$exon_species\t$exon_overlap_species\t$exon_nonoverlap_species\t$intergenic_species\t$intergenic_overlap_species\t$intergenic_nonoverlap_species\t$noannot_species\t" >> $1.annot_table.species

echo -e "Annotation\tLTR\tnon_LTR\tDNA\trepeats\tsensemRNA\tsensemRNA_overlap_transposons\tsensemRNA_nonoverlap_transposons\tantisensemRNA\tantisensemRNA_overlap_transposons\tantisensemRNA_nonoverlap_transposons\tintron\tintron_overlap_transposons\tintron_nonoverlap_transposons\texon\texon_overlap_transposons\texon_nonoverlap_transposons\tintergenic\tintergenic_overlap_transposons\tintergenic_nonoverlap_transposons\tno annotation" > $1.annot_table.reads

echo -e "All\t$LTR_reads\t$non_LTR_reads\t$DNA_reads\t$repeat_reads\t$sense_mRNA_reads\t$sense_mRNA_overlap_reads\t$sense_mRNA_nonoverlap_reads\t$antisense_mRNA_reads\t$antisense_mRNA_overlap_reads\t$antisense_mRNA_nonoverlap_reads\t$intron_reads\t$intron_overlap_reads\t$intron_nonoverlap_reads\t$exon_reads\t$exon_overlap_reads\t$exon_nonoverlap_reads\t$intergenic_reads\t$intergenic_overlap_reads\t$intergenic_nonoverlap_reads\t$noannot_reads\t" >> $1.annot_table.reads


#--------------------------piRNA annotation-------------------------------
LTR_species=`lenrangeselector $1.xkxh.transposon.mapper2.LTR.uniq.reads 23 29 | wc -l`
LTR_reads=`lenrangeselector $1.xkxh.transposon.mapper2.LTR.uniq.reads 23 29 | sumcol+ 2`

non_LTR_species=`lenrangeselector $1.xkxh.transposon.mapper2.non_LTR.uniq.reads 23 29 | wc -l`
non_LTR_reads=`lenrangeselector $1.xkxh.transposon.mapper2.non_LTR.uniq.reads 23 29 | sumcol+ 2`

DNA_species=`lenrangeselector $1.xkxh.transposon.mapper2.DNA.uniq.reads 23 29 | wc -l`
DNA_reads=`lenrangeselector $1.xkxh.transposon.mapper2.DNA.uniq.reads 23 29 | sumcol+ 2`

repeat_species=`lenrangeselector $1.xkxh.transposon.mapper2.repeats.uniq.reads 23 29 | wc -l`
repeat_reads=`lenrangeselector $1.xkxh.transposon.mapper2.repeats.uniq.reads 23 29 | sumcol+ 2`

sense_mRNA_species=`lenrangeselector $1.xkxh.snmRNA.mapper.uniq.reads 23 29 |wc -l`
sense_mRNA_reads=`lenrangeselector $1.xkxh.snmRNA.mapper.uniq.reads 23 29 | sumcol+ 2`

sense_mRNA_overlap_species=`lenrangeselector $1.xkxh.snmRNA.overlaptrans.uniq.reads 23 29 |wc -l`
sense_mRNA_overlap_reads=`lenrangeselector $1.xkxh.snmRNA.overlaptrans.uniq.reads 23 29 | sumcol+ 2`

let "sense_mRNA_nonoverlap_species=sense_mRNA_species-sense_mRNA_overlap_species"
let "sense_mRNA_nonoverlap_reads=sense_mRNA_reads-sense_mRNA_overlap_reads"

antisense_mRNA_species=`lenrangeselector $1.xkxh.antimRNA.mapper.uniq.reads 23 29 |wc -l`
antisense_mRNA_reads=`lenrangeselector $1.xkxh.antimRNA.mapper.uniq.reads 23 29 | sumcol+ 2`

antisense_mRNA_overlap_species=`lenrangeselector $1.xkxh.antimRNA.overlaptrans.uniq.reads 23 29 |wc -l`
antisense_mRNA_overlap_reads=`lenrangeselector $1.xkxh.antimRNA.overlaptrans.uniq.reads 23 29 | sumcol+ 2`

let "antisense_mRNA_nonoverlap_species=antisense_mRNA_species-antisense_mRNA_overlap_species"
let "antisense_mRNA_nonoverlap_reads=antisense_mRNA_reads-antisense_mRNA_overlap_reads"

intron_species=`lenrangeselector $1.xkxh.match2_all.out.intron.uniq.reads 23 29 |wc -l`
intron_reads=`lenrangeselector $1.xkxh.match2_all.out.intron.uniq.reads 23 29 | sumcol+ 2`

intron_overlap_species=`lenrangeselector $1.xkxh.match2_all.out.intron.overlaptrans.uniq.reads 23 29 |wc -l`
intron_overlap_reads=`lenrangeselector $1.xkxh.match2_all.out.intron.overlaptrans.uniq.reads 23 29 | sumcol+ 2`

let "intron_nonoverlap_species=intron_species-intron_overlap_species"
let "intron_nonoverlap_reads=intron_reads-intron_overlap_reads"

exon_species=`lenrangeselector $1.xkxh.match2_all.out.exon.uniq.reads 23 29 |wc -l`
exon_reads=`lenrangeselector $1.xkxh.match2_all.out.exon.uniq.reads 23 29 | sumcol+ 2`

exon_overlap_species=`lenrangeselector $1.xkxh.match2_all.out.exon.overlaptrans.uniq.reads 23 29 |wc -l`
exon_overlap_reads=`lenrangeselector $1.xkxh.match2_all.out.exon.overlaptrans.uniq.reads 23 29 | sumcol+ 2`

let "exon_nonoverlap_species=exon_species-exon_overlap_species"
let "exon_nonoverlap_reads=exon_reads-exon_overlap_reads"

intergenic_species=`lenrangeselector $1.xkxh.match2_all.out.intergenic.uniq.reads 23 29 |wc -l`
intergenic_reads=`lenrangeselector $1.xkxh.match2_all.out.intergenic.uniq.reads 23 29 | sumcol+ 2`

intergenic_overlap_species=`lenrangeselector $1.xkxh.match2_all.out.intergenic.overlaptrans.uniq.reads 23 29 |wc -l`
intergenic_overlap_reads=`lenrangeselector $1.xkxh.match2_all.out.intergenic.overlaptrans.uniq.reads 23 29 | sumcol+ 2`

let "intergenic_nonoverlap_species=intergenic_species-intergenic_overlap_species"
let "intergenic_nonoverlap_reads=intergenic_reads-intergenic_overlap_reads"

noannot_species=`lenrangeselector $1.xkxh.match2_all.out.noannot.uniq.reads 23 29 |wc -l`
noannot_reads=`lenrangeselector $1.xkxh.match2_all.out.noannot.uniq.reads 23 29 | sumcol+ 2`

echo -e "PiRNA\t$LTR_species\t$non_LTR_species\t$DNA_species\t$repeat_species\t$sense_mRNA_species\t$sense_mRNA_overlap_species\t$sense_mRNA_nonoverlap_species\t$antisense_mRNA_species\t$antisense_mRNA_overlap_species\t$antisense_mRNA_nonoverlap_species\t$intron_species\t$intron_overlap_species\t$intron_nonoverlap_species\t$exon_species\t$exon_overlap_species\t$exon_nonoverlap_species\t$intergenic_species\t$intergenic_overlap_species\t$intergenic_nonoverlap_species\t$noannot_species\t" >> $1.annot_table.species


echo -e "PiRNA\t$LTR_reads\t$non_LTR_reads\t$DNA_reads\t$repeat_reads\t$sense_mRNA_reads\t$sense_mRNA_overlap_reads\t$sense_mRNA_nonoverlap_reads\t$antisense_mRNA_reads\t$antisense_mRNA_overlap_reads\t$antisense_mRNA_nonoverlap_reads\t$intron_reads\t$intron_overlap_reads\t$intron_nonoverlap_reads\t$exon_reads\t$exon_overlap_reads\t$exon_nonoverlap_reads\t$intergenic_reads\t$intergenic_overlap_reads\t$intergenic_nonoverlap_reads\t$noannot_reads\t" >> $1.annot_table.reads

#-----------------------siRNA annot---------------------------
LTR_species=`lenselector $1.xkxh.transposon.mapper2.LTR.uniq.reads 21 | wc -l`
LTR_reads=`lenselector $1.xkxh.transposon.mapper2.LTR.uniq.reads 21 | sumcol+ 2`

non_LTR_species=`lenselector $1.xkxh.transposon.mapper2.non_LTR.uniq.reads 21 | wc -l`
non_LTR_reads=`lenselector $1.xkxh.transposon.mapper2.non_LTR.uniq.reads 21 | sumcol+ 2`

DNA_species=`lenselector $1.xkxh.transposon.mapper2.DNA.uniq.reads 21 | wc -l`
DNA_reads=`lenselector $1.xkxh.transposon.mapper2.DNA.uniq.reads 21 | sumcol+ 2`

repeat_species=`lenselector $1.xkxh.transposon.mapper2.repeats.uniq.reads 21 | wc -l`
repeat_reads=`lenselector $1.xkxh.transposon.mapper2.repeats.uniq.reads 21 | sumcol+ 2`

sense_mRNA_species=`lenselector $1.xkxh.snmRNA.mapper.uniq.reads 21 |wc -l`
sense_mRNA_reads=`lenselector $1.xkxh.snmRNA.mapper.uniq.reads 21 | sumcol+ 2`

sense_mRNA_overlap_species=`lenselector $1.xkxh.snmRNA.overlaptrans.uniq.reads 21 |wc -l`
sense_mRNA_overlap_reads=`lenselector $1.xkxh.snmRNA.overlaptrans.uniq.reads 21 | sumcol+ 2`

let "sense_mRNA_nonoverlap_species=sense_mRNA_species-sense_mRNA_overlap_species"
let "sense_mRNA_nonoverlap_reads=sense_mRNA_reads-sense_mRNA_overlap_reads"

antisense_mRNA_species=`lenselector $1.xkxh.antimRNA.mapper.uniq.reads 21 |wc -l`
antisense_mRNA_reads=`lenselector $1.xkxh.antimRNA.mapper.uniq.reads 21 | sumcol+ 2`

antisense_mRNA_overlap_species=`lenselector $1.xkxh.antimRNA.overlaptrans.uniq.reads 21 |wc -l`
antisense_mRNA_overlap_reads=`lenselector $1.xkxh.antimRNA.overlaptrans.uniq.reads 21 | sumcol+ 2`

let "antisense_mRNA_nonoverlap_species=antisense_mRNA_species-antisense_mRNA_overlap_species"
let "antisense_mRNA_nonoverlap_reads=antisense_mRNA_reads-antisense_mRNA_overlap_reads"

intron_species=`lenselector $1.xkxh.match2_all.out.intron.uniq.reads 21 |wc -l`
intron_reads=`lenselector $1.xkxh.match2_all.out.intron.uniq.reads 21 | sumcol+ 2`

intron_overlap_species=`lenselector $1.xkxh.match2_all.out.intron.overlaptrans.uniq.reads 21 |wc -l`
intron_overlap_reads=`lenselector $1.xkxh.match2_all.out.intron.overlaptrans.uniq.reads 21 | sumcol+ 2`

let "intron_nonoverlap_species=intron_species-intron_overlap_species"
let "intron_nonoverlap_reads=intron_reads-intron_overlap_reads"

exon_species=`lenselector $1.xkxh.match2_all.out.exon.uniq.reads 21 |wc -l`
exon_reads=`lenselector $1.xkxh.match2_all.out.exon.uniq.reads 21 | sumcol+ 2`

exon_overlap_species=`lenselector $1.xkxh.match2_all.out.exon.overlaptrans.uniq.reads 21 |wc -l`
exon_overlap_reads=`lenselector $1.xkxh.match2_all.out.exon.overlaptrans.uniq.reads 21 | sumcol+ 2`

let "exon_nonoverlap_species=exon_species-exon_overlap_species"
let "exon_nonoverlap_reads=exon_reads-exon_overlap_reads"

intergenic_species=`lenselector $1.xkxh.match2_all.out.intergenic.uniq.reads 21 |wc -l`
intergenic_reads=`lenselector $1.xkxh.match2_all.out.intergenic.uniq.reads 21 | sumcol+ 2`

intergenic_overlap_species=`lenselector $1.xkxh.match2_all.out.intergenic.overlaptrans.uniq.reads 21 |wc -l`
intergenic_overlap_reads=`lenselector $1.xkxh.match2_all.out.intergenic.overlaptrans.uniq.reads 21 | sumcol+ 2`

let "intergenic_nonoverlap_species=intergenic_species-intergenic_overlap_species"
let "intergenic_nonoverlap_reads=intergenic_reads-intergenic_overlap_reads"

noannot_species=`lenselector $1.xkxh.match2_all.out.noannot.uniq.reads 21 |wc -l`
noannot_reads=`lenselector $1.xkxh.match2_all.out.noannot.uniq.reads 21 | sumcol+ 2`

echo -e "siRNA\t$LTR_species\t$non_LTR_species\t$DNA_species\t$repeat_species\t$sense_mRNA_species\t$sense_mRNA_overlap_species\t$sense_mRNA_nonoverlap_species\t$antisense_mRNA_species\t$antisense_mRNA_overlap_species\t$antisense_mRNA_nonoverlap_species\t$intron_species\t$intron_overlap_species\t$intron_nonoverlap_species\t$exon_species\t$exon_overlap_species\t$exon_nonoverlap_species\t$intergenic_species\t$intergenic_overlap_species\t$intergenic_nonoverlap_species\t$noannot_species\t" >> $1.annot_table.species


echo -e "siRNA\t$LTR_reads\t$non_LTR_reads\t$DNA_reads\t$repeat_reads\t$sense_mRNA_reads\t$sense_mRNA_overlap_reads\t$sense_mRNA_nonoverlap_reads\t$antisense_mRNA_reads\t$antisense_mRNA_overlap_reads\t$antisense_mRNA_nonoverlap_reads\t$intron_reads\t$intron_overlap_reads\t$intron_nonoverlap_reads\t$exon_reads\t$exon_overlap_reads\t$exon_nonoverlap_reads\t$intergenic_reads\t$intergenic_overlap_reads\t$intergenic_nonoverlap_reads\t$noannot_reads\t" >> $1.annot_table.reads


mv $1.annot_table.reads $2
mv $1.annot_table.species $2
mv $1_stats_table_species $2
mv $1_stats_table_reads $2

rm $1.xkxh.transposon.mapper2.21nt*
rm $1.xkxh.transposon.mapper2.23-29nt*
rm $1.match2_all.out
rm $1.transposon.mapper
rm $1.transposon.mapper2
rm $1.transposon.mapper2.uniq.reads
rm $1.xkxh.transposon.mapper2.DNA.uniq.reads
rm $1.xkxh.transposon.mapper2.Foldback.uniq.reads
rm $1.xkxh.transposon.mapper2.IR.uniq.reads
rm $1.xkxh.transposon.mapper2.LTR.uniq.reads
rm $1.xkxh.transposon.mapper2.non_LTR.uniq.reads
rm $1.xkxh.transposon.mapper2.nonLTR.uniq.reads
rm $1.xkxh.transposon.mapper2.repeats.uniq.reads
rm $1.xkxh.transposon.mapper2.SINE.uniq.reads