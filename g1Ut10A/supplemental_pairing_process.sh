#!/bin/bash -x

#compare the g11-g23 to t11-t23 bit by bit, if equal then 1, otherwise 0
export PIPELINE_DIRECTORY=/home/wangw1/git/smallRNA_analysis
for i in `seq 1 10`; 
do
	${PIPELINE_DIRECTORY}/g1Ut10A/compair_f1_f3 AubIPuniq_w1_unox_.Ago3IPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.r${i} > AubIPuniq_w1_unox_.Ago3IPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.sum.r${i} && \
	rm -rf AubIPuniq_w1_unox_.Ago3IPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.r${i}
	${PIPELINE_DIRECTORY}/g1Ut10A/compair_f1_f3 Ago3IPuniq_w1_unox_.AubIPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.r${i} > Ago3IPuniq_w1_unox_.AubIPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.sum.r${i} && \
	rm -rf Ago3IPuniq_w1_unox_.AubIPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.r${i}
done