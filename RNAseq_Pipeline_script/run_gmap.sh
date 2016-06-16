#! /usr/bin/env bash

dbpath=/diag/home/netashawang/cloud/spu_gmap/Spur_3.1.LinearScaffold/
genome=Spur_3.1.LinearScaffold
tri_fa=$1

FILE=${1##*/}
PREFIX=${FILE%.fasta*}

gmap -D $dbpath -d $genome -k 15 -t 20 --format=psl -n 500 --quiet-if-excessive --suboptimal-score=0 --fails-as-input --split-output=${PREFIX}.gmap $tri_fa > ${PREFIX}.gmap.sam

