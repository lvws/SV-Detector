#!/usr/bin/env python
"""
用于合并Fusion Hunter 和 Husion Catch; 方便上传lims
格式：
id	Type	Gene	Gene.ID	AAChange	Chr.start	Chr.end	Ref	Alt	Hom.Het	ExonicFunc	rsID	AF	CopyNumber	DP(ref:alt)	NormalID	NormalAF	NormalDP(ref: alt)	Comment	Hotspot	Strand	CLNSIG	CLNREVSTAT	CLNDNINCL	CLINID	BRCAClass	CosmicID	Cosmic.Occurence	HGMD_class	Variant_effect_LOVD	InterVar_Class	1000g	1000gEAS	ExAC_ALL	ExAC_EAS	gnomAD_exome_ALL	gnomAD_exome_EAS	gnomAD_genome_ALL	gnomAD_genome_EAS	GS_ALL	SIFT.score	SIFT.pred	PolyPhen.score	PolyPhen.pred	CADD_phred	GERP + +_RS	OMIM_Inheritance	OMIM_Link	InterVar_ACMG	Chr.start_rule3	Chr.end_rule3	VCFKEY	HGVS3aa	HGVSg	BrkpFusion	BrkptType	TumorReferenceCount	TumorSplitReferenceCount	TumorVariantCount	TumorSplitVariantCount	Cosmic_Fusion_Counts	CC_Tumour_Types(Somatic)	Orientation	PIPELINEVERSION	DP1	DP2	AllSupportCount
FA219C0163-R108	SV	NTRK3	ETV6:NM_001987.5|NTRK3:NM_001012338.2	ETV6:exon5~NTRK3:exon15	12:12022903	15:88483984	.	.	.	.	.	22.63%	.	1772	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	5	1553	.	.	.	.	1772	1191	1558

在 NormalAF 列标记融合检出来自何处：
r0: 两个软件都有检出
r1: fc 检出
r2: fh 检出
"""
import sys,os

if len(sys.argv) < 3:
    print(sys.argv[0],"fusion_catcher.tsv\tfusion_hunter.tsv")
    sys.exit(1)

fc_file = sys.argv[1]
fh_file = sys.argv[2]

# 获取文件中突变信息
def getInfo(infile):
    dic = {}
    with open(infile,'r') as f:
        header = next(f)
        h_list = header.split('\t')
        tindex = h_list.index("NormalAF")
        si = h_list.index("Chr.start")
        ei = h_list.index("Chr.end")
        for i in f:
            t = i.split('\t')
            s = t[si]
            e = t[ei]
            dic[s+'-'+e] = [t,False]
    return header,tindex,dic

# 比较两个位点是否临近（默认 500 flank）
def checkPos(s1,s2,flank=500):
    c1,p1 = s1.split(':')
    c2,p2 = s2.split(':')
    if c1 != c2:
        return False
    if abs(int(p1) - int(p2)) < flank:
        return True
    return False


header,tindex,fc_dic = getInfo(fc_file)
header,tindex,fh_dic = getInfo(fh_file)
print(header,end='')
for ss in fc_dic:
    tag = 'r1'
    cs1,cs2 = ss.split('-')
    for hh in fh_dic:
        hs1,hs2 = hh.split('-')
        if (checkPos(cs1,hs1) and checkPos(cs2,hs2)) or (checkPos(cs1,hs2) and checkPos(cs2,hs1)):
            tag = 'r0'
            fh_dic[hh][1] = True
            break
    fc_dic[ss][0][tindex] = tag
    print('\t'.join(fc_dic[ss][0]),end='')

for hh in fh_dic:
    tag = 'r2'
    if fh_dic[hh][1]:
        continue
    hs1,hs2 = hh.split('-')
    for ss in fc_dic:
        cs1,cs2 = ss.split('-')
        if (checkPos(cs1,hs1) and checkPos(cs2,hs2)) or (checkPos(cs1,hs2) and checkPos(cs2,hs1)):
            tag = 'r0'
            break
    fh_dic[hh][0][tindex] = tag
    print('\t'.join(fh_dic[hh][0]),end='')




