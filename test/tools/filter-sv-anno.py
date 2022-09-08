#!/usr/bin/env python
from ctypes import alignment
import os,sys
from collections import defaultdict
import pysam

"""
# annot.tsv
chr1	pos1	str1	chr2	pos2	str2	gene1	transcript1	site1	kinase_domain1	gene2	transcript2	site2	kinase_domain2	fusion	Cosmic_Fusion_Counts	repName-repClass-repFamily:-site1	repName-repClass-repFamily:-site2	CC_Chr_Band	CC_Tumour_Types(Somatic)	CC_Cancer_Syndrome	CC_Mutation_Type	CC_Translocation_Partner	DGv_Name-DGv_VarType-site1	DGv_Name-DGv_VarType-site2

# panel bed
1	7721786	7721845	WWTR1_C4_CAMTA1_C8
1	7723413	7723472	WWTR1_C4_CAMTA1_C9
1	7723566	7723639	WWTR1_C3a_CAMTA1_C9,WWTR1_C2_I2b_CAMTA1_C9,WWTR1_C2_I2a_CAMTA1_C9
"""

if len(sys.argv) < 4:
	print(sys.argv[0],"annot.txt\tpanel.bed\traw-sv.tsv\tbam")
	sys.exit(1)

anno_file = sys.argv[1]
panel_file = sys.argv[2]
sv_file = sys.argv[3]
bam_file = sys.argv[4] # 用于获取融合点位的depth
Alignment = pysam.AlignmentFile(bam_file,'rb')

sample = os.path.basename(sv_file).split('.')[0]
# get bed regin infos(assume low->high)
bed_dic = defaultdict(lambda:[])


# 获取融合点位的depth
def depthGet(site):
	chrom,pos = site.split(':')
	return Alignment.count(chrom,int(pos),int(pos)+1)

with open(panel_file,'r') as f:
	for i in f:
		t = i.strip().split('\t')
		bed_dic[t[0]].append([int(t[1]),int(t[2])])

#check fusion pos is in bed region
def checkRegin(chrom,pos,bed_dic,flank=200):
	if bed_dic == {}:
		return True
	pos_int = int(pos)
	for s,e in bed_dic[chrom]:
		if pos_int > e + flank:
			continue
		elif pos_int + flank > s:
			return True
	return False 

fheader = "id,Type,Gene,Gene.ID,AAChange,Chr.start,Chr.end,Ref,Alt,Hom.Het,ExonicFunc,rsID,AF,CopyNumber,DP(ref:alt),NormalID,NormalAF,NormalDP(ref: alt),Comment,Hotspot,Strand,CLNSIG,CLNREVSTAT,CLNDNINCL,CLINID,BRCAClass,CosmicID,Cosmic.Occurence,HGMD_class,Variant_effect_LOVD,InterVar_Class,1000g,1000gEAS,ExAC_ALL,ExAC_EAS,gnomAD_exome_ALL,gnomAD_exome_EAS,gnomAD_genome_ALL,gnomAD_genome_EAS,GS_ALL,SIFT.score,SIFT.pred,PolyPhen.score,PolyPhen.pred,CADD_phred,GERP + +_RS,OMIM_Inheritance,OMIM_Link,InterVar_ACMG,Chr.start_rule3,Chr.end_rule3,VCFKEY,HGVS3aa,HGVSg,BrkpFusion,BrkptType,TumorReferenceCount,TumorSplitReferenceCount,TumorVariantCount,TumorSplitVariantCount,Cosmic_Fusion_Counts,CC_Tumour_Types(Somatic),Orientation,PIPELINEVERSION,DP1,DP2,AllSupportCount".split(',')
print('\t'.join(fheader))

# 获取原始信息
"""
break_pos1	break_pos2	soft_map_dr(p1);soft_map_dr(p2)	gene_p1_dr;gene_p2_dr	support_p1_dr;support_p2_dr	gene_p1;gene_p2	split_p1:split_p2:discordant reads	p1_cigars;p2_cigars;total_reads
9:131456290	3:155705836	+;+	+;+	u;d	SET:NM_003011:e7:25;NA	31:41:0	20:16	AGATGATGATGATGATGAA:GAGGAGGAAGGATTAGAAGATATTGACGAAGAAGGGGATGAGGATGAAGGTGAAGAAGA
"""
dic = {}
with open(sv_file,'r') as f:
	lines = next(f).strip().split('\t')
	ig1 = lines.index('break_pos1')
	ig2 = lines.index('break_pos2')
	imdr = lines.index('soft_map_dr(p1);soft_map_dr(p2)')
	isp = lines.index('split_p1:split_p2:discordant reads')
	ispl = lines.index('p1_cigars;p2_cigars;total_reads')
	for i in f:
		t = i.strip().split('\t')
		if len(t) < 8:
			break
		g1,g2,mdr = t[:3]
		span = t[isp].split(':')[-1]
		spl = t[ispl]
		seq = t[ispl+3].replace(':','*')
		key = g1 + '-' + g2
		dic[key] = [mdr,span,spl,seq]

# 拼凑AAChange
"""
可能的形式：
Exon 2 of TPST1(+)
Intron of JAZF1(-):19Kb after exon 2
5'-UTR of ADGRL3(+):296Kb before coding start
IGR: 73Kb before FAM92B(-)
IGR: 199Kb before CCND2(+)

转变的形式：
IGR (downstream GRB14)~ALK:exon20
RET:exon112~IGR (upstream TRIB1)
AKT1:exon11~ROS1:exon35
ROS1:intron34~TRIB1:exon10
EZR:exon9~ROS1:exon34
FGFR3:exon18~TACC3:exon11
FGFR1:5'UTR~IGR (upstream TRMT12)
"""
def makeAAChange(gene,site):
	slist = site.split(' ')
	if slist[0] == "Exon":
		return gene + ':exon' + slist[1]
	elif slist[0] == "Intron":
		num = '1'
		if slist[-1].isdigit():
			num = slist[-1]
		return gene + ':intron' + num
	elif slist[0] == 'IGR:':
		tag = '(downstream ' + gene + ')'
		if slist[-1].endswith('(+)'):
			tag = '(upstream ' + gene + ')'
		return "IGR " + tag
	else:
		return gene + ':' + slist[0]


"""
id	Type	Gene	Gene.ID	AAChange	Chr.start	Chr.end	Ref	Alt	Hom.Het	ExonicFunc	rsID	AF	CopyNumber	DP(ref:alt)	NormalID	NormalAF	NormalDP(ref: alt)	Comment	Hotspot	Strand	CLNSIG	CLNREVSTAT	CLNDNINCL	CLINID	BRCAClass	CosmicID	Cosmic.Occurence	HGMD_class	Variant_effect_LOVD	InterVar_Class	1000g	1000gEAS	ExAC_ALL	ExAC_EAS	gnomAD_exome_ALL	gnomAD_exome_EAS	gnomAD_genome_ALL	gnomAD_genome_EAS	GS_ALL	SIFT.score	SIFT.pred	PolyPhen.score	PolyPhen.pred	CADD_phred	GERP + +_RS	OMIM_Inheritance	OMIM_Link	InterVar_ACMG	Chr.start_rule3	Chr.end_rule3	VCFKEY	HGVS3aa	HGVSg	BrkpFusion	BrkptType	TumorReferenceCount	TumorSplitReferenceCount	TumorVariantCount	TumorSplitVariantCount	Cosmic_Fusion_Counts	CC_Tumour_Types(Somatic)	Orientation	PIPELINEVERSION	DP1	DP2	AllSupportCount
FA219C0163-R108	SV	NTRK3	ETV6:NM_001987.5|NTRK3:NM_001012338.2	ETV6:exon5~NTRK3:exon15	12:12022903	15:88483984	.	.	.	.	.	22.63%	.	1772	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	5	1553	.	.	.	.	1772	1191	1558
"""
# 过滤panel 内的融合，并修改，表示方式。
with open(anno_file,'r') as f:
	header = next(f)
	#print(header,end="")
	header_lst = header.split('\t')
	for i in f:
		t = i.split('\t')
		if t[header_lst.index('gene1')] == t[header_lst.index('gene2')]:
			continue
		chr1 = t[header_lst.index('chr1')]
		pos1 = t[header_lst.index('pos1')]
		chr2 = t[header_lst.index('chr2')]
		pos2 = t[header_lst.index('pos2')]
		if checkRegin(chr1,pos1,bed_dic) or checkRegin(chr2,pos2,bed_dic):
			info_dic = {}
			g1 = chr1 + ':' + pos1
			g2 = chr2 + ':' + pos2
			mdr,span,spl,seq = dic[g1+'-'+g2]
			seq1,seq2 = seq.split('*')
			seqlen = max(len(seq1),len(seq2))
			spll,splr = spl.split(':')
			spc = int(spll) + int(splr)
			asc = int(span) + spc
			md1,md2 = mdr.split(';')
			gene1 = t[header_lst.index('gene1')]
			gene2 = t[header_lst.index('gene2')]
			trans1 = t[header_lst.index('transcript1')]
			trans2 = t[header_lst.index('transcript2')]
			site1 = t[header_lst.index('site1')]
			site2 = t[header_lst.index('site2')]

			# 获取depth 信息，计算丰度
			dp1 = str(depthGet(g1))
			dp2 = str(depthGet(g2))
			dp = max(int(dp1),int(dp2))
			if dp1 > dp2:
				af = str(round(float(spll)/int(dp1)*100,2)) + "%"
			else:
				af = str(round(float(splr)/int(dp2)*100,2)) + "%"

			info_dic['id'] = sample
			info_dic['Type'] = 'SV'
			info_dic['Gene'] = gene1
			info_dic['Gene.ID'] = "{gene1}:{trans1}|{gene2}:{trans2}".format(**locals())
			info_dic['AAChange'] = makeAAChange(gene1,site1) + '~' + makeAAChange(gene2,site2)
			info_dic['Chr.start'] = g1
			info_dic['Chr.end'] = g2
			info_dic['AF'] = af
			info_dic['DP(ref:alt)'] = str(dp)
			info_dic['TumorVariantCount'] = str(span)
			info_dic['TumorSplitVariantCount'] = str(spc)
			info_dic['DP1'] = dp1
			info_dic['DP2'] = dp2
			info_dic['AllSupportCount'] = str(asc)

			info_list = []
			for h in fheader:
				info_list.append(info_dic.get(h,'.'))
			
			print('\t'.join(info_list))

		else:
			sys.stderr.write("Not in bed regin: " + i)

Alignment.close()
