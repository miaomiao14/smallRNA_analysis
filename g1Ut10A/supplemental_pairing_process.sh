#!/bin/bash -x
for i in `seq 1 10`; 
do
	./compair_f1_f3 AubIPuniq_w1_unox_.Ago3IPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.r${i} > AubIPuniq_w1_unox_.Ago3IPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.sum.r${i} && \
	rm -rf AubIPuniq_w1_unox_.Ago3IPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.r${i}
	./compair_f1_f3 Ago3IPuniq_w1_unox_.AubIPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.r${i} > Ago3IPuniq_w1_unox_.AubIPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.sum.r${i} && \
	rm -rf Ago3IPuniq_w1_unox_.AubIPuniq_w1_unox_.10.prefix.UA_VA.ppseq.txt.shuffled.r${i}
done