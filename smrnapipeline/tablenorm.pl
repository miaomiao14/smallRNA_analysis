#!/usr/bin/perl
$insert=$ARGV[0];

open IN, "$insert.norm.bed";
while(<IN>) {
chomp; split(/\t/);
$key="$_[0]:$_[1]-$_[2]($_[3])";
$NTM{$key}=$_[6];
}

open IN, "$insert.intergenic.mapper";
while(<IN>) {
chomp; split(/\t/);
$intergenic{$_[2]}++;
}

open IN, "$insert.intron.mapper";
while(<IN>) {
chomp; split(/\t/);
$intron{$_[2]}++;
}
open IN, "$insert.exon.mapper";
while(<IN>) {
chomp; split(/\t/);
$exon{$_[2]}++;
}
open IN, "$insert.xkxh.transposon.mapper2";
while(<IN>) {
chomp; split(/\t/);
$transposon{$_[2]}++;
}


open IN, "$insert.xkxh.match2_all.out";
while(<IN>) {
chomp; split(/\t/);
$key="$_[2]:$_[4]-$_[5]($_[3])";
if (!exists $intergenic{$key} && !exists $intron{$key} && !exists $exon{$key} && !exists $transposon{$key}) {
$noannot_species+=1/$NTM{$key}; $noannot_reads+=$_[1]/$NTM{$key};
if (length($_[0]) >=23 && length($_[0]) <=29) { $noannot_species_piRNA+=1/$NTM{$key}; $noannot_reads_piRNA+=$_[1]/$NTM{$key};}
if (length($_[0]) ==21) { $noannot_species_siRNA+=1/$NTM{$key}; $noannot_reads_siRNA+=$_[1]/$NTM{$key};}
}
else {
$intergenic{$key}=0 if (!exists $intergenic{$key});
$intron{$key}=0 if (!exists $intron{$key});
$exon{$key}=0 if (!exists $exon{$key});
$transposon{$key}=0 if (!exists $transposon{$key});
$total=$intergenic{$key}+$intron{$key}+$exon{$key}+$transposon{$key};

if($intron{$key} !=0 && $total != 0) { $intron_species+=1/$NTM{$key}/$total*$intron{$key}; $intron_reads+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_species+=1/$NTM{$key}/$total*$exon{$key}; $exon_reads+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_species+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_reads+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_species+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_reads+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
if  (length($_[0]) >=23 && length($_[0]) <=29) {
if($intron{$key} !=0 && $total != 0) { $intron_species_piRNA+=1/$NTM{$key}/$total*$intron{$key}; $intron_reads_piRNA+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_species_piRNA+=1/$NTM{$key}/$total*$exon{$key}; $exon_reads_piRNA+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_species_piRNA+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_reads_piRNA+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_species_piRNA+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_reads_piRNA+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
}
if  (length($_[0]) ==21) {
if($intron{$key} !=0 && $total != 0) {$intron_species_siRNA+=1/$NTM{$key}/$total*$intron{$key}; $intron_reads_siRNA+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_species_siRNA+=1/$NTM{$key}/$total*$exon{$key}; $exon_reads_siRNA+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_species_siRNA+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_reads_siRNA+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_species_siRNA+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_reads_siRNA+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
}
if ($transposon{$key}) {
if($intron{$key} !=0 && $total != 0) {$intron_overlaptrans_species+=1/$NTM{$key}/$total*$intron{$key}; $intron_overlaptrans_reads+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_overlaptrans_species+=1/$NTM{$key}/$total*$exon{$key}; $exon_overlaptrans_reads+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_overlaptrans_species+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_overlaptrans_reads+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_overlaptrans_species+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_overlaptrans_reads+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
if  (length($_[0]) >=23 && length($_[0]) <=29) {
if($intron{$key} !=0 && $total != 0) {$intron_overlaptrans_species_piRNA+=1/$NTM{$key}/$total*$intron{$key}; $intron_overlaptrans_reads_piRNA+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_overlaptrans_species_piRNA+=1/$NTM{$key}/$total*$exon{$key}; $exon_overlaptrans_reads_piRNA+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_overlaptrans_species_piRNA+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_overlaptrans_reads_piRNA+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_overlaptrans_species_piRNA+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_overlaptrans_reads_piRNA+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
}
if  (length($_[0]) ==21) {
if($intron{$key} !=0 && $total != 0) {$intron_overlaptrans_species_siRNA+=1/$NTM{$key}/$total*$intron{$key}; $intron_overlaptrans_reads_siRNA+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_overlaptrans_species_siRNA+=1/$NTM{$key}/$total*$exon{$key}; $exon_overlaptrans_reads_siRNA+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_overlaptrans_species_siRNA+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_overlaptrans_reads_siRNA+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_overlaptrans_species_siRNA+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_overlaptrans_reads_siRNA+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
}
}
else{
if($intron{$key} !=0 && $total != 0) {$intron_nonoverlaptrans_species+=1/$NTM{$key}/$total*$intron{$key}; $intron_nonoverlaptrans_reads+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_nonoverlaptrans_species+=1/$NTM{$key}/$total*$exon{$key}; $exon_nonoverlaptrans_reads+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_nonoverlaptrans_species+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_nonoverlaptrans_reads+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_nonoverlaptrans_species+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_nonoverlaptrans_reads+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
if  (length($_[0]) >=23 && length($_[0]) <=29) {
if($intron{$key} !=0 && $total != 0) {$intron_nonoverlaptrans_species_piRNA+=1/$NTM{$key}/$total*$intron{$key}; $intron_nonoverlaptrans_reads_piRNA+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_nonoverlaptrans_species_piRNA+=1/$NTM{$key}/$total*$exon{$key}; $exon_nonoverlaptrans_reads_piRNA+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_nonoverlaptrans_species_piRNA+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_nonoverlaptrans_reads_piRNA+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_nonoverlaptrans_species_piRNA+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_nonoverlaptrans_reads_piRNA+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
}
if  (length($_[0]) ==21) {
if($intron{$key} !=0 && $total != 0) { $intron_nonoverlaptrans_species_siRNA+=1/$NTM{$key}/$total*$intron{$key}; $intron_nonoverlaptrans_reads_siRNA+=$_[1]/$NTM{$key}/$total*$intron{$key};}
if($exon{$key} !=0 && $total != 0) {$exon_nonoverlaptrans_species_siRNA+=1/$NTM{$key}/$total*$exon{$key}; $exon_nonoverlaptrans_reads_siRNA+=$_[1]/$NTM{$key}/$total*$exon{$key};}
if($intergenic{$key} !=0 && $total != 0) {$intergenic_nonoverlaptrans_species_siRNA+=1/$NTM{$key}/$total*$intergenic{$key}; $intergenic_nonoverlaptrans_reads_siRNA+=$_[1]/$NTM{$key}/$total*$intergenic{$key};}
if($transposon{$key} !=0 && $total != 0) {$transposon_nonoverlaptrans_species_siRNA+=1/$NTM{$key}/$total*$transposon{$key}; $transposon_nonoverlaptrans_reads_siRNA+=$_[1]/$NTM{$key}/$total*$transposon{$key};}
}
}
}
}

open OUT, ">$insert.annot_table.norm.reads";
print OUT  "Annotation\ttransposon\tintron\tintron_overlap_transposons\tintron_nonoverlap_transposons\texon\texon_overlap_transposons\texon_nonoverlap_transposons\tintergenic\tintergenic_overlap_transposons\tintergenic_nonoverlap_transposons\tno_annotation\n"; 

print OUT "All\t$transposon_reads\t$intron_reads\t$intron_overlaptrans_reads\t$intron_nonoverlaptrans_reads\t$exon_reads\t$exon_overlaptrans_reads\t$exon_nonoverlaptrans_reads\t$intergenic_reads\t$intergenic_overlaptrans_reads\t$intergenic_nonoverlaptrans_reads\t$noannot_reads\n";
print OUT "piRNA\t$transposon_reads_piRNA\t$intron_reads_piRNA\t$intron_overlaptrans_reads_piRNA\t$intron_nonoverlaptrans_reads_piRNA\t$exon_reads_piRNA\t$exon_overlaptrans_reads_piRNA\t$exon_nonoverlaptrans_reads_piRNA\t$intergenic_reads_piRNA\t$intergenic_overlaptrans_reads_piRNA\t$intergenic_nonoverlaptrans_reads_piRNA\t$noannot_reads_piRNA\n";
print OUT "siRNA\t$transposon_reads_siRNA\t$intron_reads_siRNA\t$intron_overlaptrans_reads_siRNA\t$intron_nonoverlaptrans_reads_siRNA\t$exon_reads_siRNA\t$exon_overlaptrans_reads_siRNA\t$exon_nonoverlaptrans_reads_siRNA\t$intergenic_reads_siRNA\t$intergenic_overlaptrans_reads_siRNA\t$intergenic_nonoverlaptrans_reads_siRNA\t$noannot_reads_siRNA\n";
print "CHECK:reads_total=",$transposon_reads+$intron_reads+$exon_reads+$intergenic_reads+$noannot_reads,"\n";

open OUT, ">$insert.annot_table.norm.species";
print OUT "Annotation\ttransposon\tintron\tintron_overlap_transposons\tintron_nonoverlap_transposons\texon\texon_overlap_transposons\texon_nonoverlap_transposons\tintergenic\tintergenic_overlap_transposons\tintergenic_nonoverlap_transposons\tno_annotation\n";


print OUT "All\t$transposon_species\t$intron_species\t$intron_overlaptrans_species\t$intron_nonoverlaptrans_species\t$exon_species\t$exon_overlaptrans_species\t$exon_nonoverlaptrans_species\t$intergenic_species\t$intergenic_overlaptrans_species\t$intergenic_nonoverlaptrans_species\t$noannot_species\n";
print OUT "piRNA\t$transposon_species_piRNA\t$intron_species_piRNA\t$intron_overlaptrans_species_piRNA\t$intron_nonoverlaptrans_species_piRNA\t$exon_species_piRNA\t$exon_overlaptrans_species_piRNA\t$exon_nonoverlaptrans_species_piRNA\t$intergenic_species_piRNA\t$intergenic_overlaptrans_species_piRNA\t$intergenic_nonoverlaptrans_species_piRNA\t$noannot_species_piRNA\n";
print OUT "siRNA\t$transposon_species_siRNA\t$intron_species_siRNA\t$intron_overlaptrans_species_siRNA\t$intron_nonoverlaptrans_species_siRNA\t$exon_species_siRNA\t$exon_overlaptrans_species_siRNA\t$exon_nonoverlaptrans_species_siRNA\t$intergenic_species_siRNA\t$intergenic_overlaptrans_species_siRNA\t$intergenic_nonoverlaptrans_species_siRNA\t$noannot_species_siRNA\n";

print "CHECK:species_total=",$transposon_species+$intron_species+$exon_species+$intergenic_species+$noannot_species,"\n";

`mv $insert.annot_table.norm.reads $ARGV[1]/`;
`mv $insert.annot_table.norm.species $ARGV[1]/`;