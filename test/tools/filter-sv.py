#!/usr/bin/env python
import os,sys

if len(sys.argv) < 2:
	print(sys.argv[0],"fusion_raw.tsv")
	sys.exit(1)

with open(sys.argv[1],'r') as f:
	for i in f:
		if i.startswith('break'):
			print(i,end='')
			break
	for i in f:
		t = i.strip().split('\t')
		if len(t) < 3:
			break
		gene_infos = t[5].split(';')
		gene1 = gene_infos[0].split(':')[0]
		gene2 = gene_infos[1].split(':')[0]
		extron1_off = 0
		extron2_off = 0
		if gene1 != "NA":
			extron1_off = int(gene_infos[0].split(':')[-1])
		if gene2 != "NA":
			extron2_off = int(gene_infos[1].split(':')[-1])
		ref_b1,ref_b2,dis = map(lambda x : int(x), t[6].split(':'))
		ref_b1_u ,ref_b2_u  = map(lambda x:int(x),t[7].split(':'))
		if(gene1 == gene2 and gene1 != "NA"):
			sys.stderr.write("same gene !\t" + i)
			continue
		if (gene1 == "NA" or gene2 == "NA" or max(extron1_off,extron2_off) > 10):
			if(ref_b1_u + ref_b2_u >= 10 and min(ref_b1_u,ref_b2_u) >= 3):
				print(i,end='')
			else:
				sys.stderr.write("out extron split read low!{ref_b1_u}:{ref_b2_u}\t".format(**locals())+i)
		elif(dis >= 6 or min(ref_b1_u,ref_b2_u) >= 2) :
			print(i,end='')
		else:
			sys.stderr.write("in extron split read low!{dis};{ref_b1_u}:{ref_b2_u}\t".format(**locals())+i)

