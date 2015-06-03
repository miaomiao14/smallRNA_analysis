#!/bin/sh
FA=$1
TYPE=$2

CPU=${3:-8}
OUTDIR=$4
prefix=$5
#FILE=${1##*/}
#prefix=`basename $FILE _best_candidates.eclipsed_orfs_removed.pep`

PFAM=/diag/home/netashawang/indexes/databases/Pfam-A.hmm
DB=/diag/home/netashawang/indexes/databases/uniprot_sprot.fasta

#get fasta from bed
#gtf_to_fasta $GTF $GENOMEFA ${GTF}.fa

[ ! -d $OUTDIR ] && mkdir -p $OUTDIR 
if [[ $TYPE == "fasta" ]]
then 
#rename the fasta sequence from fastafrombed
#if use -name then this step is unnecessary
#awk 'BEGIN{OFS="\t";}{if(/>/){print ">"$3,$2;} else {print $0;}}' $FA > $OUTDIR/${prefix}.fa
#run transdecoder
ln -s $FA $OUTDIR/${prefix}.fa
/diag/home/netashawang/src/trinityrnaseq_r2013-02-25/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl -t $OUTDIR/${prefix}.fa -S && \
mv $OUTDIR/best_candidates.eclipsed_orfs_removed.pep $OUTDIR/${prefix}_best_candidates.eclipsed_orfs_removed.pep && \
PEP=$OUTDIR/${prefix}_best_candidates.eclipsed_orfs_removed.pep
elif [[ $TYPE == "pep" ]]
then
FILE=${1##*/}
[ ! -f ${OUTDIR}/$FILE ] && ln -s $FA ${OUTDIR}/${FILE}
PEP=${OUTDIR}/${FILE}
fi


[ ! -f ${OUTDIR}/${prefix}_TrinotateBlast.out  ] && ~/src/ncbi-blast-2.2.28+/bin/blastp -query $PEP -db $DB -num_threads $CPU -max_target_seqs 1 -outfmt 6 -out ${OUTDIR}/${prefix}_TrinotateBlast.out
[ ! -f ${OUTDIR}/${prefix}_TrinotatePFAM.out  ] && hmmscan --cpu $CPU --domtblout ${OUTDIR}/${prefix}_TrinotatePFAM.out1 $PFAM $PEP > ${OUTDIR}/pfam.log
grep -v "#" ${OUTDIR}/${prefix}_TrinotatePFAM.out1 >${OUTDIR}/${prefix}_TrinotatePFAM.out
[ ! -f ${OUTDIR}/${prefix}_signalp.out ] && signalp -f short -n ${OUTDIR}/${prefix}_signalp.out $PEP
[ ! -f ${OUTDIR}/${prefix}_tmhmm.out ] && tmhmm --short < $PEP > ${OUTDIR}/${prefix}_tmhmm.out

[ ! -f ${OUTDIR}/TrinityFunctional.db ] && ln -s /diag/home/netashawang/indexes/databases/TrinityFunctional.db ${OUTDIR}/TrinityFunctional.db
/diag/home/netashawang/src/trinityrnaseq_r2013-02-25/Analysis/Trinotate/Trinotate.pl LOAD_transdecoder ${OUTDIR}/$FILE && \
/diag/home/netashawang/src/trinityrnaseq_r2013-02-25/Analysis/Trinotate/Trinotate.pl LOAD_blast ${OUTDIR}/${prefix}_TrinotateBlast.out && \
/diag/home/netashawang/src/trinityrnaseq_r2013-02-25/Analysis/Trinotate/Trinotate.pl LOAD_pfam ${OUTDIR}/${prefix}_TrinotatePFAM.out && \
/diag/home/netashawang/src/trinityrnaseq_r2013-02-25/Analysis/Trinotate/Trinotate.pl LOAD_tmhmm ${OUTDIR}/${prefix}_tmhmm.out && \
/diag/home/netashawang/src/trinityrnaseq_r2013-02-25/Analysis/Trinotate/Trinotate.pl LOAD_signalp ${OUTDIR}/${prefix}_signalp.out && \
/diag/home/netashawang/src/trinityrnaseq_r2013-02-25/Analysis/Trinotate/Trinotate.pl report > ${OUTDIR}/${prefix}_trinotate_annotation_report.xls
