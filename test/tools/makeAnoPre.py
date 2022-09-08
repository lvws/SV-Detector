#!/usr/bin/env python
import os,sys

if len(sys.argv) < 2:
	print(sys.argv[0],"filter.tsv")
	sys.exit(1)

print("chr1\tpos1\tstr1\tchr2\tpos2\tstr2")
with open(sys.argv[1],'r') as f:
	next(f)
	for i in f:
		t = i.strip().split('\t')
		c1,p1 = t[0].split(':')
		c2,p2 = t[1].split(':')
		s1,s2 = t[3].split(';')
		ss1 = '1'
		ss2 = '1'
		if s1 == '-':
			ss1 = '0'
		if s2 == '-':
			ss2 = '0'
		print('\t'.join([c1,p1,ss1,c2,p2,ss2]))
		
