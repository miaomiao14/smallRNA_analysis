#!/bin/sh
dir=$1
fast=$2
paraFile=${RANDOM}${RANDOM}.para
for i in `ls $dir/*.fastq`
do
prefix=`basename $i .fastq`
#echo $prefix
if ( -z $fast)
{
	echo -e "awk -v left=\"${prefix}_left.fq\" -v right=\"${prefix}_right.fq\" '{\
	print \$0\"/1\" >>left; \
	print \$0\"/2\" >>right;\
	getline; \
	L=length(\$0); \
	print substr(\$0,0,rshift(L,1)) >> left; \
	print substr(\$0,rshift(L,1)+1) >> right; \
	getline; \
	getline; \
	print \"+\\\n\"substr(\$0,0,rshift(L,1)) >> left; \
	print \"+\\\n\"substr(\$0,rshift(L,1)+1) >> right; \
	}' $i" >> $paraFile
}
else
	{
	echo -e "awk 'NR%2==1 { print \$0\"/1\" } ; NR%2==0 { print substr(\$0,1,length(\$0)/2) }' $i >${prefix}_left.fq" >> $paraFile
	echo -e "awk 'NR%2==1 { print \$0\"/2\" } ; NR%2==0 { print substr(\$0,length(\$0)/2+1) }' $i >${prefix}_right.fq" >> $paraFile
}
fi
done
CPU=`wc -l $paraFile`
ParaFly -c $paraFile -CPU $CPU