
##input files is ppseq.txt from pp8_v24_v2.pl
input=$1

outdir=$3
##anchor point, extend 100 nt on both end

##trans
grep trans $input |cut -f7 \
	| awk 'BEGIN{OFS="\t"}{split($1,a,","); if(a[3]=="+"){start=a[2]+10;strand="-";left=start-100;right=start+100; print a[1],left,right,"*","*",strand} else {start=a[2]-9;strand="+";print a[1],left,right,"*","*",strand}}' \
		|sort -u >${outdir}/${input}.trans.bed
grep "+" ${input}.trans.bed >${outdir}/${input}.trans.plus.bed
grep "-" ${input}.trans.bed >${outdir}/${input}.trans.minus.bed
##cis
grep cis $input |cut -f7 \
	| awk 'BEGIN{OFS="\t"}{split($1,a,","); if(a[3]=="+"){start=a[2]+10;strand="-"; left=start-100;right=start+100; print a[1],left,right,"*","*",strand} else {start=a[2]-9;strand="+";print a[1],left,right,"*","*",strand}}' \
		|sort -u >${outdir}/${input}.cis.bed &
grep "+" ${input}.cis.bed >${outdir}/${input}.trans.plus.bed
grep "-" ${input}.cis.bed >${outdir}/${input}.trans.minus.bed

##use the above bed intervals as chrsize

normbed=$2
normFac=1
zcat $normbed |sed 1d | awk -F "\t" -v norm=$normFac '{OFS="\t"; if($4=="+" && length($5)>=23) print $1,$2-1,$3,$6","$7","norm/1}' >${outdir}/${normbed%.norm.bed.gz}.bed.plus.pirna 
zcat $normbed |sed 1d | awk -F "\t" -v norm=$normFac '{OFS="\t"; if($4=="-" && length($5)>=23) print $1,$2-1,$3,$6","$7","norm/1}' >${outdir}/${normbed%.norm.bed.gz}.bed.plus.pirna 


###bedItemOverlapCountWithScore does not work on zlab3
bedtools sort -i ${outdir}/${normbed%.norm.bed.gz}.bed.plus.pirna | bedItemOverlapCountWithScore dm3 stdin -chromSize=${outdir}/${input}.trans.plus.bed  >${outdir}/${normbed%.norm.bed.gz}.bed.plus.sorted.pirna 
bedtools sort -i ${outdir}/${normbed%.norm.bed.gz}.bed.minus.pirna | bedItemOverlapCountWithScore dm3 stdin -chromSize=${outdir}/${input}.trans.minus.bed >${outdir}/${normbed%.norm.bed.gz}.bed.minus.sorted.pirna 


bedtools sort -i ${outdir}/${normbed%.norm.bed.gz}.bed.plus.pirna | bedItemOverlapCountWithScore dm3 stdin -chromSize=${outdir}/${input}.cis.plus.bed  >${outdir}/${normbed%.norm.bed.gz}.bed.cis.sorted.pirna 

