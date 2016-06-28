#!/bin/bash

STAR --runMode genomeGenerate --genomeDir /home/sgeadmin/data/mf5.0/refseq/ --genomeFastaFiles /home/sgeadmin/data/mf5.0/refseq/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fna --sjdbGTFfile/home/sgeadmin/data/mf5.0/refseq/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.gtf  --runThreadN 8 --limitGenomeGenerateRAM 640000000000
