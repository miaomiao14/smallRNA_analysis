for i in VOY*; do
        echo -ne "$i\n";
		
		##use array in awk
		awk 'BEGIN{OFS="\t"}{if ($1~/guide/) { a+=$4/$5;b[$2]+=$4/$5;} else { c+=$4/$5; e[$2]+=$4/$5;}} END { n=asorti(b,bd); for (ii = 1; ii <= n; ii++) {r=b[bd[ii]]/a; if(bd[ii]>=0) {print "guide","N+"bd[ii],r;}else {print "guide","N"bd[ii],r;} } m=asorti(e,ed); for (j = 1; j <= m; j++) {s=e[ed[j]]/c; if(ed[j]>=0){print "passenger","N+"ed[j],s;}else{print "passenger","N"ed[j],s;} }}' $i/final.bed2 >>percentage.txt;
		awk -v sample=$i 'BEGIN{OFS="\t";}{if ($1~/guide/) { a+=$4/$5;b[$2]+=$4/$5;} else { c+=$4/$5; e[$2]+=$4/$5;}} END {ORS="\t"; {print sample;} n=asorti(b,bd); for (ii = 1; ii <= n; ii++) {r=b[bd[ii]]/a; if(bd[ii]==1) {print "N+"bd[ii]":"r;}if(bd[ii]==0) {print "N:"r;}if(bd[ii]==-1) {print "N"bd[ii]":"r;} } ORS="";print "\n"} ' $i/final.bed2 >>percentage.subset.txt;
done
