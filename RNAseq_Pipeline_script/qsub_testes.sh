qsub -P diag -pe make 24 -q highmem.q -V -b y -N "trinity_SeaUrchin_Testes" -cwd /diag/home/pkulance/softwares/trinityrnaseq_r2012-10-05/Trinity.pl --seqType fq --JM 100G --left /diag/home/pkulance/scratch/RNASeq/sea_urchin/130116_I282_FCD1LCPACXX_L7_CHKPEI12120266_1.fq --right /diag/home/pkulance/scratch/RNASeq/sea_urchin/130116_I282_FCD1LCPACXX_L7_CHKPEI12120266_2.fq --SS_lib_type RF --output /diag/home/netashawang/scratch/trinity_testes_spu --CPU 24 --min_kmer_cov 1 --inchworm_cpu 24 --bflyHeapSpaceInit 1G --bflyHeapSpaceMax 800G --bflyCPU 20 --bflyGCThreads 4
