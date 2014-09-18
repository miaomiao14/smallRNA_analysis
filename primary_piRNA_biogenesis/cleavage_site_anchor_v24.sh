
##input files is ppseq.txt from pp8_v24_v2.pl
#use prefix 16 instead of prefix 10
input=$1
outdir=$2

##use target coordinates as anchor point
##cleavage site as anchor point,  
##trans
grep trans $input | awk 'BEGIN{OFS="\t"}{split($7,a,","); if(a[3]=="+"){start=a[2]+10;strand="-"; print a[1],start-1,start,$9,$5*$9,strand} else {start=a[2]-10;strand="+";print a[1],start,start+1,$9,$5*$9,strand}}' > ${outdir}/${input}.trans &&
sort -k1,1 -k2,2 -k3,3 -k6,6 ${outdir}/${input}.trans |/home/wangw1/wangwj2bin/bedtools groupby -i stdin -g 1,2,3,6 -c 4,5 -o sum,sum >${outdir}/${input}.trans.groupby && 
sort -k 6,6rn ${outdir}/${input}.trans.groupby > ${outdir}/${input}.trans.bed &
#rm ${outdir}/${input}.trans ${outdir}/${input}.trans.groupby &
##cis
grep cis $input | awk 'BEGIN{OFS="\t"}{split($7,a,","); if(a[3]=="+"){start=a[2]+10;strand="-";  print a[1],start-1,start,$9,$5*$9,strand} else {start=a[2]-10;strand="+";print a[1],start,start+1,$9,$5*$9,strand}}' | gsort -k1,1 -k2,2 -k3,3 -k6,6 -T $PWD --parallel=8 | /home/wangw1/wangwj2bin/bedtools groupby -i stdin -g 1,2,3,6 -c 4,5 -o sum,sum | gsort -k 6,6rn -T $PWD --parallel=8 >${outdir}/${input}.cis.bed &

##can also try to use guide coordinates(but the same piwi protein) as anchor point



