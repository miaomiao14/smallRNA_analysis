#!/usr/bin/env bash
for i in `ls *.VA.pp`
do
UA_VA_zscore.pl $i
done