#! /usr/bin/env bash
adaptor_trimmer_phred_dumper $1 TGGAATTC 0 18 5 | awk -v trimlen=3 '{print $0; getline; print substr($1,1,length($1)-trimlen); getline; print "+"; getline; print substr($1,1,length($1)-trimlen);}' | gzip > ${1%fq}trimmed.fq.gz
gzip $1
