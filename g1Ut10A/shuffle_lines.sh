#!/bin/bash -x

#shuffle the lines in a file

for i in `seq 1 10`; do
	awk 'BEGIN{getline;}{if($1=="trans") print $4"\t"$5}' $1 | shuf > ${1}.targetpiRNAs.shuffled
	awk 'BEGIN{getline;}{if($1=="trans") print $8"\t"$9}' $1 >  ${1}.f89
	paste ${1}.targetpiRNAs.shuffled ${1}.f89 > ${1}.shuffled.r${i}
	rm -rf ${1}.targetpiRNAs.shuffled ${1}.f89
done