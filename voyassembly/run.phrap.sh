for i in 1 2 3; do
    cat cyno01-P3-SCt-0${i}-Seq1.seq cyno01-P3-SCt-0${i}-Seq2.seq >  cyno01-P3-SCt-0${i}-Seq12.seq
    ../phrap  cyno01-P3-SCt-0${i}-Seq12.seq
    cat cyno01-P3-SCt-0${i}-Seq12.seq.contigs cyno01-P3-SCt-0${i}-Seq5.seq >  cyno01-P3-SCt-0${i}-Seq125.seq
    ../phrap cyno01-P3-SCt-0${i}-Seq125.seq
done

