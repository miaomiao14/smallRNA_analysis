#!/bin/bash


OUTDIR=/home/wangw1/src/circos-0.56/fly/etc
[ -f ${OUTDIR}/file.conf ] && rm ${OUTDIR}/file.conf
STATFILE=CIRCOSBIN.${RANDOM}.stat
for file in "$@"
do
CIRCOSBIN=$file
filename=${CIRCOSBIN##*/}
awk -v f=$filename 'BEGIN{OFS="\t"}{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4< min) {min=$4}; total+=$4; count+=1;} END {print f,total/count, min, max}' $CIRCOSBIN >>${OUTDIR}/${STATFILE}
done

maxv=`awk 'BEGIN{OFS="\t"}{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4< min) {min=$4};}END {print min}' ${OUTDIR}/${STATFILE}`
maxvalue=`awk -v m=$maxv 'BEGIN{print m/2}'`
echo "
#karyotype = ../data/karyotype.drosophila.hires.dm3.txt
karyotype = ../data/karyotype.fly.all.dm3.txt
chromosomes_units           = 1000000
chromosomes_display_default = yes
#chromosomes_scale = chrx:2
#chromosomes_order = chrx,chr2l,chr3l,chr2r,chr3r,chr4
#relative order
#chromosomes_order = chr1,chr4,‐,chr3,chr5
#chromosomes = chrx:8-12;chr2l:4-11;chr3R:9-21
#chromosomes_breaks = -chr2l:11-19
	<plots>
" >>${OUTDIR}/file.conf
num_vars=$#
count=0
radius_inc=$(awk -v a=${num_vars} 'BEGIN{print 0.6/a}')
base=0.4
for CIRCOSBIN in "$@"
do
filename=${CIRCOSBIN##*/}
rr0=$(awk -v a=$base -v b=$count -v c=$radius_inc 'BEGIN{print a+b*c}')
rr1=$(awk -v a=$base -v b=$count -v c=$radius_inc 'BEGIN{print a+(b+1)*c}')
if [[ $filename =~ "minus" ]] || [[ $filename =~ "antisense" ]]
then 

echo "
	<plot>
	type  = line
	min   = 0
	max   = $maxvalue
	orientation = in
	file  = ${CIRCOSBIN}
	r0    = ${rr0}r
	r1    = ${rr1}r
	color = vvdred
	thickness = 0.75
	</plot> " >>${OUTDIR}/file.conf
fi

if [[ $filename =~ "plus" ]] || [[ $filename =~ "\.sense" ]]
then

echo "

	<plot>
	type  = line
	min   = 0
	max   = $maxvalue
	file  = ${CIRCOSBIN}
	r0    = ${rr0}r
	r1    = ${rr1}r
	color = vvdblue
	thickness = 0.75
	</plot> " >>${OUTDIR}/file.conf
fi


count=$(($count+1))

done

echo "
	</plots>
<<include ideogram.conf>>
<<include ticks.conf>>

<image>
<<include ../../etc/image.conf>>
</image>
<<include housekeeping.conf>> " >>${OUTDIR}/file.conf