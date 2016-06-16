#!/usr/bin/env python

import sys

fh = open(sys.argv[1])
while True:
	line = fh.readline()
	if not line: break
	if line[0] == '>':
		name = line[1:].rstrip()
		print name + ".fa"
		out = open(name + ".fa", 'w')
	else:
		out.write(">" + name + "\n" + line + '\n')
	
