for i in *Seq1.seq; do
    	name=${i/cyno01-P3-Cau-/};

	n=${name/-Seq1.seq/};	
	#echo $n;
    cat cyno01-P3-Cau-${n}-Seq1.seq cyno01-P3-Cau-${n}-Seq2.seq >  cyno01-P3-Cau-${n}-Seq12.seq
    phrap  cyno01-P3-Cau-${n}-Seq12.seq
    cat cyno01-P3-Cau-${n}-Seq12.seq.contigs cyno01-P3-Cau-${n}-Seq5.seq >  cyno01-P3-Cau-${n}-Seq125.seq
    phrap cyno01-P3-Cau-${n}-Seq125.seq
done

