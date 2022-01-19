#!/usr/bin/env python
import os,sys,pickle
from pyfaidx import Faidx
from collections import defaultdict

if len(sys.argv) < 4:
    print(sys.argv[0],"reference(hs37d5.fa)\tseed length\tindexName")
    sys.exit(1)

fasta = sys.argv[1]
seedLen = int(sys.argv[2])
indexName = sys.argv[3]

refer = Faidx(fasta)
chrom_list = ["1"]
# chrom_list = list("123456789") + ["10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
dic = {}

with open(fasta+".fai",'r') as f:
    for i in f:
        t = i.strip().split('\t')
        chrom = t[0]
        length = int(t[1])
        if chrom in chrom_list:
            print("Procesing ",chrom)
            n = 1
            while n + seedLen - 1 < length:
                seq = refer.fetch(chrom,n,n+seedLen-1).seq
                if not 'N' in seq:
                    if seq in dic:
                        dic[seq].append(chrom+':'+str(n))
                    else:
                        dic[seq] = [chrom+':'+str(n)]
                n += 1



with open(indexName,'wb') as f:
    pickle.dump(dic,f)

