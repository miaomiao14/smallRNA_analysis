#!/usr/bin/env python
import sys
from collections import defaultdict
def count(file):
    vec = defaultdict(float)
    for line in open(file):
        tokens = line.rstrip().split()
        bit = int(tokens[0])
        for i in range(12,-1,-1):
            if bit & 1<<i:
                # print "{0:b}\t{1}".format(bit, tokens[1])
                vec[12-i] += float(tokens[1])
            else:
                break
    for i in vec:
        print "{0}\t{1}".format(i, vec[i])
if __name__=="__main__":
    count(sys.argv[1])