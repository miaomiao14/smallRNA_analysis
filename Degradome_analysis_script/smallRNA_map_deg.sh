#! /bin/bash -x

# get genomic mappers of deg
# warning: do NOT use bedtools bamtofastq, which messes  up \1 and \2
for i in *bam; do submitsge 8 "$i" "samtools bam2fq $i | awk '{a=substr(\$1,length(\$1)); head=\$1;getline;seq=\$1; getline;getline; qual=\$1; if (a==1) {printf \"%s\\n%s\\n+\\n%s\\n\",head,seq,qual}}' > ${i%bam}r1.fq"; done
for i in *bam; do submitsge 8 "$i" "samtools bam2fq $i | awk '{a=substr(\$1,length(\$1)); head=\$1;getline;seq=\$1; getline;getline; qual=\$1; if (a==2) {printf \"%s\\n%s\\n+\\n%s\\n\",head,seq,qual}}' > ${i%bam}r2.fq"; done


for i in *r1.fq; do submitsge 2 $i "catFastqUnique $i ${i%%_*} > ${i%fq}fa && bowtie-build ${i%fq}fa ${i%%_*}"; done
for i in *r2.fq; do submitsge 2 $i "catFastqUnique $i ${i%%_*}ctrl > ${i%fq}fa && bowtie-build ${i%fq}fa ${i%%_*}ctrl"; done

#
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.ago3_CD.x_rRNA.dm3.sorted.F0x100.r1
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.ago3_CD.x_rRNA.dm3.sorted.F0x100.r2
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.ago3_mut.x_rRNA.dm3.sorted.F0x100.r1
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.ago3_mut.x_rRNA.dm3.sorted.F0x100.r2
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.ago3_WT.x_rRNA.dm3.sorted.F0x100.r1
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.ago3_WT.x_rRNA.dm3.sorted.F0x100.r2
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.aub_CD.x_rRNA.dm3.sorted.F0x100.r1
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.aub_CD.x_rRNA.dm3.sorted.F0x100.r2
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.aub_mut.x_rRNA.dm3.sorted.F0x100.r1
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.aub_mut.x_rRNA.dm3.sorted.F0x100.r2
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.aub_WT.x_rRNA.dm3.sorted.F0x100.r1
/home/hanb/scratch/CD/degradomeBam/Phil.Deg.aub_WT.x_rRNA.dm3.sorted.F0x100.r2	

/home/hanb/scratch/cd/Phil.Deg.ago3_CD/degIndex/Phil.Deg.ago3_CD.r1		Phil.aubvasAgo3CDrescue.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.ago3_mut/degIndex/Phil.Deg.ago3_mut.r1	Phil.Ago3MutsWW.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.ago3_WT/degIndex/Phil.Deg.ago3_WT.r1		Phil.aubvasAgo3WTrescue.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.aub_CD/degIndex/Phil.Deg.aub_CD.r1		Phil.AubCDrescue.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.aub_mut/degIndex/Phil.Deg.aub_mut.r1		Phil.AubMutsWW.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.aub_WT/degIndex/Phil.Deg.aub_WT.r1		Phil.AubWTrescue.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.ago3_CD/degIndex/Phil.Deg.ago3_CD.r2		Phil.aubvasAgo3CDrescue.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.ago3_mut/degIndex/Phil.Deg.ago3_mut.r2	Phil.Ago3MutsWW.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.ago3_WT/degIndex/Phil.Deg.ago3_WT.r2		Phil.aubvasAgo3WTrescue.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.aub_CD/degIndex/Phil.Deg.aub_CD.r2		Phil.AubCDrescue.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.aub_mut/degIndex/Phil.Deg.aub_mut.r2		Phil.AubMutsWW.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa
/home/hanb/scratch/cd/Phil.Deg.aub_WT/degIndex/Phil.Deg.aub_WT.r2		Phil.AubWTrescue.ox.ovary.trimmed.x_rRNA.x_Hairpin.insert.x_dmel_virus.dm3v0a.al.fa




# map
for i in `cat indexList`; do submitsge 8 "${i##*/}" "bowtie -f -a --best --strata -v 0 -p 8 -S $i ~/scratch/smallRNA/fasta/RF.w1118.ovary.unox.trim.dm3a.unique.l23mer.fa | samtools view -bS - > ${i##*/}.w1118.bam"; done

# convert bamt to bed
for i in `ls | grep bam`; do echo "bedtools bamtobed -i $i > ${i%bam}bed" >> bamToBed.para; done
ParaFly -c bamToBed.para -CPU 20

# convert bed to bed2
for i in *.bed; do echo "0.convert.sh $i" >> 0.para; done
ParaFly -c 0.para -CPU 20

# get + and - count
for i in *bed2; do echo "1.plus_minus.sh $i" >> 1.para; done
ParaFly -c 1.para -CPU 20

for i in *bed2; do echo "4.countAll.sh  $i" >> 4.para; done
ParaFly -c 4.para -CPU 20

# counting sense mapper 5' end
awk '{if ($6=="+") {ct[$2%100]+=1.0/$5}}END{for (a=0;a<=99;++a){printf "%d\t%.2f\n",a,ct[a]}}' $1 > ${1}.overlapCount.sense5Prime
# counting sense mapper 3' end
awk '{if ($6=="+") {ct[$3%100]+=1.0/$5}}END{for (a=0;a<=99;++a){printf "%d\t%.2f\n",a,ct[a]}}' $1 > ${1}.overlapCount.sense3Prime
# counting antisense mapper 5' end
awk '{if ($6=="-") {ct[$3%100]+=1.0/$5}}END{for (a=0;a<=99;++a){printf "%d\t%.2f\n",a,ct[a]}}' $1 > ${1}.overlapCount.antisense5Prime
# counting antisense mapper 3' end
awk '{if ($6=="-") {ct[$2%100]+=1.0/$5}}END{for (a=0;a<=99;++a){printf "%d\t%.2f\n",a,ct[a]}}' $1 > ${1}.overlapCount.antisense3Prime
# draw

for i in `ls | grep overlapCount | grep -v sh`; do echo "~/nearline/small_RNA_Pipeline/bin/Rscript  ~/nearline/small_RNA_Pipeline/bin/draw_pp.R $i ${i}" >> drawPP.para; done
ParaFly -c drawPP.para -CPU 20



PREFIX=w1118.species.s.5p
PDFs=`ls | grep pdf | grep sense5Prime`
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PREFIX}.pdf ${PDFs} && rm -rf $PDFs

PREFIX=w1118.species.s.3p
PDFs=`ls | grep pdf | grep sense3Prime`
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PREFIX}.pdf ${PDFs} && rm -rf $PDFs

PREFIX=w1118.species.as.5p
PDFs=`ls | grep pdf | grep antisense5Prime`
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PREFIX}.pdf ${PDFs} && rm -rf $PDFs

PREFIX=w1118.species.as.3p
PDFs=`ls | grep pdf | grep antisense3Prime`
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PREFIX}.pdf ${PDFs} && rm -rf $PDFs


# change 2.overlap.sh 
rm -rf drawPP.para.completed 2.para.completed 
ParaFly -c 2.para -CPU 20
ParaFly -c drawPP.para -CPU 20
PREFIX=w1118.species.-.5p

for i in *.bed; do echo "0.convert.sh $i" >> 0.para; done
ParaFly -c 0.para -CPU 20
for i in *bed2; do echo "4.countAll.sh  $i" >> 4.para; done
ParaFly -c 4.para -CPU 20
for i in `ls | grep overlapCount | grep -v sh`; do echo "~/nearline/small_RNA_Pipeline/bin/Rscript  ~/nearline/small_RNA_Pipeline/bin/draw_pp.R $i ${i}" >> drawPP.para; done
ParaFly -c drawPP.para -CPU 20
PREFIX=w1118.species.s.5p
PDFs=`ls | grep pdf | grep sense5Prime`
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PREFIX}.pdf ${PDFs} && rm -rf $PDFs

PREFIX=w1118.species.s.3p
PDFs=`ls | grep pdf | grep sense3Prime`
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PREFIX}.pdf ${PDFs} && rm -rf $PDFs

PREFIX=match.species.s.5p
PDFs=`ls | grep pdf | grep sense5Prime`
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PREFIX}.pdf ${PDFs} && rm -rf $PDFs

PREFIX=match.species.s.3p
PDFs=`ls | grep pdf | grep sense3Prime`
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${PREFIX}.pdf ${PDFs} && rm -rf $PDFs







