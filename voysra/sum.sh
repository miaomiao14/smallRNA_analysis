for i in VOY*; do
	echo -ne "$i\t";
	awk '{a+=$2}END{printf "%.1f\t", a}' $i/*trim.insert;
	awk '{printf "%.1f\t", $1}' $i/*hsa.hairpin.count;
	awk '{printf "%.1f\t", $1}' $i/*x_hsa.hairpin*VOY*.count;
	awk '{if ($1~/guide/) { a+=$4/$5; if ($2==0) { b+=$4/$5; if ($3==0) c+=$4/$5 } } \
		else {d+=$4/$5; if ($2==0) {e+=$4/$5; if ($3==0) f+=$4/$5}}}END{printf "%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n", a,b,c,d,e,f}' $i/final.bed2;
done
