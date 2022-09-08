#!/usr/bin/env python
import os,sys

if len(sys.argv) < 2:
	print(sys.argv[0],"annot.tsv")
	sys.exit(1)

with open(sys.argv[1],'r') as f:
	print(next(f),end='')
	for i in f:
		t = i.strip().split('\t')
		if t[6] == t[10]:
			continue
		print(i,end='')
