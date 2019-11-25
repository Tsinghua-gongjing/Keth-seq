import gj
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
import sys, os
from nested_dict import nested_dict
import pandas as pd
import numpy as np
from pyfasta import Fasta
import os
import re
import normRT_corr

def read_RG4_regions(RG4_txt=None):
	RG4_txt = RG4_txt if RG4_txt is not None else '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Guo_2016_Science_RG4/mESC_RG4_full_table1_6140RG4.ensembl.txt'
	RG4_regions_dict = nested_dict(1, list)
	with open(RG4_txt, 'r') as TXT:
		for line in TXT:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			arr = line.split('\t')
			""" # mESC_RG4.txt
			ensembl_tx_id = arr[10]
			start = int(arr[11])
			end = int(arr[12])
			"""
			ensembl_tx_id = arr[11]
			start = int(arr[12])
			end = int(arr[13])
			RG4_regions_dict[ensembl_tx_id].append([start, end])
	print RG4_regions_dict['ENSMUST00000079496']
	return RG4_regions_dict.to_dict()

def read_icshape_out(out=None, pureID=1):
	gj.printFuncRun('read_icshape_out')
	gj.printFuncArgs()
	out_dict = nested_dict()
	with open(out, 'r') as OUT:
		for line in OUT:
			line = line.strip()
			if not line or line.startswith('#'): continue
			arr = line.split('\t')
			tx_id = arr[0]
			if pureID:
				tx_id = tx_id.split('.')[0]
			length = int(arr[1])
			rpkm = float(arr[2]) if arr[2] != '*' else arr[2]
			reactivity_ls = arr[3:]
			out_dict[tx_id]['tx_id'] = tx_id
			out_dict[tx_id]['length'] = length
			out_dict[tx_id]['rpkm'] = rpkm
			out_dict[tx_id]['reactivity_ls'] = reactivity_ls
	gj.printFuncRun('read_icshape_out')
	return out_dict

def read_fa(fa='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'):
	gj.printFuncRun('read_fa')
	gj.printFuncArgs()
	fa_dict1 = Fasta(fa, key_fn=lambda key:key.split("\t")[0])
	fa_dict = {i.split()[0].split('.')[0]:j[0:] for i,j in fa_dict1.items()}
	print fa_dict.keys()[0:3]
	gj.printFuncRun('read_fa')
	return fa_dict

def icshape_in_RG4(out=None, write_noInRG4_fa=1, keep_G_only=0):
	out = out if out is not None else '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.T0t0.out'
	out_dict = read_icshape_out(out)
	fa_dict = read_fa('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/hg38/transcriptome/hg38_transcriptome.fa')
	RG4_regions_dict = read_RG4_regions('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Kit_2016_NatureMethod_RG4/nmeth.3965-S4-tableS3.convert.ensembl.txt')

	stat = out.replace('.out', '.inRG4.stat.txt')
	savefn = out.replace('.out', '.inRG4.txt')
	SAVEFN = open(savefn, 'w')

	stat_dict = nested_dict(1, int)
	stat_dict['condition'] = out.split('/')[-1].split('_')[2]
	stat_dict['T_t_cutoff'] = out.split('/')[-1].split('.')[1]

	if write_noInRG4_fa:
		savefn_fa = out.replace('.out', '.notinRG4.fa')
		SAVEFN_FA = open(savefn_fa, 'w')

	for tx in out_dict:
		stat_dict['tx'] += 1
		if RG4_regions_dict.has_key(tx):
			stat_dict['in_RG4_tx'] += 1
			RG4_len = sum([int(j)-int(i) for i,j in RG4_regions_dict[tx]])
			stat_dict['in_RG4_tx_LenRG4'] += RG4_len

			for start,end in RG4_regions_dict[tx]:
				# print >>SAVEFN, '\t'.join(map(str, [tx, start, end, ','.join(out_dict[tx]['reactivity_ls'][start:end])]))
				stat_dict['in_RG4_tx_region'] += 1
				if keep_G_only:
					G_reactivity_ls = [i for i,j in zip(out_dict[tx]['reactivity_ls'][start:end], fa_dict[tx][start:end]) if j.upper() == 'G']
					G_null_percent = len([i for i in G_reactivity_ls if i == 'NULL']) / float(len(out_dict[tx]['reactivity_ls'][start:end]))
					if G_null_percent >= 0.6:
						region_gini = -1
					else:
						region_gini = gj.gini(G_reactivity_ls, mode='mean_reactivity', null_pct=1)
						out_dict[tx]['reactivity_ls'][start:end] = [i if j.upper() == 'G' else 'NULL' for i,j in zip(out_dict[tx]['reactivity_ls'][start:end], fa_dict[tx][start:end])]
				else:
					region_gini = gj.gini(out_dict[tx]['reactivity_ls'][start:end], mode='gini', null_pct=0.6)
				if region_gini > 0:
					stat_dict['in_RG4_tx_region_validStructure'] += 1
					print >>SAVEFN, '\t'.join(map(str, [tx, start, end, region_gini, ','.join(out_dict[tx]['reactivity_ls'][start:end]), fa_dict[tx][start:end]]))
		else:
			stat_dict['not_in_RG4_tx'] += 1
			if write_noInRG4_fa:
				print >>SAVEFN_FA, '>%s'%(tx)
				print >>SAVEFN_FA, fa_dict[tx][0:]
	SAVEFN.close()
	if write_noInRG4_fa:
		SAVEFN_FA.close()

	print stat_dict
	stat_df = pd.DataFrame.from_dict(stat_dict, orient='index')
	stat_df.columns = [out.split('/')[-1]]
	print stat_df.T

	return stat_df.T


def seq_base_ratio(seq=None, base='G'):
	base_ratio = seq.count(base) / float(len(seq))
	return base_ratio

def read_fa2(fa=None):
	fa = fa if fa is not None else '/Share/home/zhangqf7/gongjing/zebrafish/script/zhangting/paris_RBP/Figure-RBP/combine/combine_utr3.fa'
	fa_dict = {}
	with open(fa, 'r') as FA:
		for n,line in enumerate(FA):
			if not line or line.startswith('#'): continue
			if line.startswith('>'):
				seq_id = line.strip().replace('>', '')
				fa_dict[seq_id] = FA.next().strip()
	# print fa_dict['NM_001001847:1134-1184']
	return fa_dict



def RG4_compare(RG4_1=None, RG4_2=None, savefn_prefix=None):
	# T_t_cutoff = 'T1t200'
	# if RG4_1 is None:
	# 	RG4_1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_mRNA_kethoxal.%s.inRG4.txt'%(T_t_cutoff)
	# if RG4_2 is None:
	# 	RG4_2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4_kethoxal.%s.inRG4.txt'%(T_t_cutoff)
	df1 = pd.read_csv(RG4_1, header=None, sep='\t')
	df2 = pd.read_csv(RG4_2, header=None, sep='\t')
	header_ls1 = ['tx', 'start', 'end', 'mRNA gini', 'reactivity_ls', 'seq']
	header_ls2 = ['tx', 'start', 'end', 'G4 gini', 'reactivity_ls', 'seq']
	df1.columns = header_ls1
	df2.columns = header_ls2
	df1['region'] = ['%s:%s-%s'%(t,s,e) for t,s,e in zip(df1['tx'],df1['start'],df1['end'])]
	df2['region'] = ['%s:%s-%s'%(t,s,e) for t,s,e in zip(df2['tx'],df2['start'],df2['end'])]

	# savefn = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/19-04-03/mRNA_G4.%s.RG4.region_overlap.png'%(T_t_cutoff)
	savefn = '%s.RG4.region_overlap.png'%(savefn_prefix)
	gj.venn3plot(mode='string',subsets_ls=[set(df1['region'].value_counts().keys()), set(df2['region'].value_counts().keys())],labels_ls=['mRNA', 'G4'],title_str=None,save_fn=savefn,axis=None)

	df_join = pd.merge(df1, df2, on=['tx', 'start', 'end'], how='inner')
	df_join_more_structure = df_join[df_join['G4 gini'] > df_join['mRNA gini']]
	print "all: %s, G4 more structure: %s, less structure: %s"%(df_join.shape[0], df_join_more_structure.shape[0], df_join.shape[0]-df_join_more_structure.shape[0])

	fig,ax=plt.subplots(figsize=(6,6))
	# region_label = ['ENST00000356674:920-950']
	# region_color = ['black' if i in region_label else 'grey' for i in df_join['region_x']]
	# df_join.plot(kind='scatter', x='mRNA gini', y='G4 gini', ax=ax, color='region_color')

	df_join.plot(kind='scatter', x='mRNA gini', y='G4 gini', ax=ax, color='#8172B2')
	# locs, labels = [0.6, 0.7, 0.8, 0.9, 1.0], [0.6, 0.7, 0.8, 0.9, 1.0]
	# plt.xticks(locs, labels)
	# plt.yticks(locs, labels)
	#x = np.linspace(*ax.get_xlim())

	# plot for subregion
	# ax.plot([0.45, 1],[0.45, 1],ls='--', color='black')
	# ax.axis('square')
	# plt.xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
	# plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
	# ax.set_xlim(0.45,1)
	# ax.set_ylim(0.45,1)

	ax.plot([0, 1],[0, 1],ls='--', color='black')
	plt.xticks([0,  0.2,  0.4,  0.6,  0.8,  1.0])
	plt.yticks([0,  0.2,  0.4,  0.6,  0.8,  1.0])
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)

	plt.tight_layout()
	# savefn = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/19-04-03/mRNA_G4.%s.RG4.region_gini.scatter.pdf'%(T_t_cutoff)
	savefn = '%s.RG4.region_gini.scatter.pdf'%(savefn_prefix)
	plt.savefig(savefn)
	plt.close()

	# return
	print df_join_more_structure.columns
	print df_join_more_structure
	# plot_more_structure(df=df_join_more_structure, extend=20)
	df_join.to_csv(savefn.replace('.pdf', '.txt'), header=True, index=False, sep='\t')
	return df_join_more_structure

def plot_more_structure(df, extend=20):
	mRNA_T0t0 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.T0t0.out'
	mRNA_T0t20 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.T0t20.out'
	G4_T0t0 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4_kethoxal.T0t0.out'
	G4_T0t20 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4_kethoxal.T0t20.out'
	mRNA_T0t0_out = read_icshape_out(mRNA_T0t0)
	mRNA_T0t20_out = read_icshape_out(mRNA_T0t20)
	G4_T0t0_out = read_icshape_out(G4_T0t0)
	G4_T0t20_out = read_icshape_out(G4_T0t20)
	kethoxal_noPDS_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-K-1_S4_L006_R1_001.rt'
	kethoxal_PDS_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-K-2_S5_L006_R1_001.rt'
	kethoxal_noPDS_rt_dict = normRT_corr.read_rt(rt=kethoxal_noPDS_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)
	kethoxal_PDS_rt_dict = normRT_corr.read_rt(rt=kethoxal_PDS_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)

	fa_dict = read_fa('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/hg38/transcriptome/hg38_transcriptome.fa')

	save_dir = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-09-09_10_library_G4-glucose_time/mRNA_vs_G4'
	for tx,start,end,gini1,gini2 in zip(df['tx'], df['start'], df['end'], df['mRNA gini'], df['G4 gini']):
		seq = fa_dict[tx][start:end+extend]
		reactivity_ls1 = map(float, mRNA_T0t0_out[tx]['reactivity_ls'][start:end+extend])
		reactivity_ls2 = map(float, G4_T0t0_out[tx]['reactivity_ls'][start:end+extend])
		rt_ls1 = map(float, kethoxal_noPDS_rt_dict[tx]['RTstop'][start:end+extend])
		rt_ls2 = map(float, kethoxal_PDS_rt_dict[tx]['RTstop'][start:end+extend])
		savefn = save_dir+'/'+'%s.%s-%s.png'%(tx, start, end)
		title = '%s:%s-%s, %.3f vs %.3f'%(tx, start, end, gini1, gini2)
		ylim1 = max(reactivity_ls1+reactivity_ls2)
		ylim2 = max(rt_ls1+rt_ls2)
		plot_multiple_shape_track(seq=seq, reactivity_ls_ls=[reactivity_ls1, reactivity_ls2, rt_ls1, rt_ls2], highlight_base_ls=[], title=title, savefn=savefn, ylim_ls=[ylim1, ylim1, ylim2, ylim2])

def plot_multiple_shape_track(seq=None, reactivity_ls_ls=None, highlight_base_ls=None, title='', savefn=None, ylim_ls=None, base='ATCGU', set_non_G_0=1):
    if seq is None:
        seq = 'ATCGGTCGAA' * 8
    if reactivity_ls_ls is None:
        reactivity_ls_ls = [[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] * 8,
                            [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1] * 8]
    if highlight_base_ls is None:
        highlight_base_ls = [4, 5, 6]
    if savefn is None:
        savefn = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_region_candidate/5UTR/test.shape.png'
    # if ylim_ls is None:
    	# ylim_ls = [6,6,6]
    reactivity_ls_ls = [ [j for j,b in zip(i, list(seq)) if b in base] for i in reactivity_ls_ls]
    seq = ''.join([b for b in list(seq) if b in base])
    color_ls = ['#C44E52' if n +
                1 in highlight_base_ls else '#4C72B0' for n, i in enumerate(list(seq))]
    fig, ax = plt.subplots(len(reactivity_ls_ls), 1,
                           sharex=True, sharey=False, figsize=(14, 12))
    for n, reactivity_ls in enumerate(reactivity_ls_ls):
        reactivity_ls = [-0.1 if i ==
                         'NULL' else float(i) for i in reactivity_ls]
        if set_non_G_0:
        	reactivity_ls = [i if j == 'G' else -0.1 for i,j in zip(reactivity_ls, list(seq))]
        	print reactivity_ls
        ax[n].bar(xrange(len(reactivity_ls)),
                  reactivity_ls, width=1.0, color=color_ls)
        if ylim_ls is not None and ylim_ls[n]:
        	ax[n].set_ylim(0,ylim_ls[n])
        plt.subplots_adjust(hspace=.001)
    plt.xticks(range(0, len(reactivity_ls)), list(seq))
    plt.suptitle(title)
    plt.subplots_adjust(wspace=0.02)
    plt.tight_layout()
    # savefn = '/Share/home/zhangqf7/gongjing/zebrafish/result/dynamic_region_candidate/5UTR/test.shape.png'
    plt.savefig(savefn)
    plt.close()
    return fig

def in_vivo_in_vitro_compare(txt1=None, txt2=None):
	if txt1 is None: txt1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-09-09_10_library_G4-glucose_time/mRNA_G4.T0t20.RG4.region_gini.scatter.txt'
	if txt2 is None: txt2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/19-07-05/in_vitro_G4_mRNA.T1t200.RG4.region_gini.scatter.txt'
	df1 = pd.read_csv(txt1, header=0, sep='\t')
	df2 = pd.read_csv(txt2, header=0, sep='\t')
	df1 = df1[df1['G4 gini'] > df1['mRNA gini']]
	df2 = df2[df2['G4 gini'] > df2['mRNA gini']]
	df1['region'] = ['%s:%s-%s'%(t,s,e) for t,s,e in zip(df1['tx'], df1['start'], df1['end'])]
	df2['region'] = ['%s:%s-%s'%(t,s,e) for t,s,e in zip(df2['tx'], df2['start'], df2['end'])]
	savefn = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/19-07-05/in_vivo_in_vitro_compare.morestructure.overlap.png'
	gj.venn3plot(mode='string',subsets_ls=[set(df1['region'].value_counts().keys()), set(df2['region'].value_counts().keys())],labels_ls=['in vivo', 'in vitro'],title_str=None,save_fn=savefn,axis=None)

	df_join = pd.merge(df1, df2, on=['tx', 'start', 'end'], how='inner')
	df_join.to_csv(savefn.replace('.png', '.txt'), header=True, index=False, sep='\t')

def in_vivo_in_vitro_compare_reactivity(extend=20):
	txt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/19-07-05/in_vivo_in_vitro_compare.morestructure.overlap.txt'
	df = pd.read_csv(txt, header=0, sep='\t')
	in_vivo_G4_T0t20 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4_kethoxal.T0t20.out'
	in_vivo_mRNA_T0t20 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.T0t20.out'
	in_vitro_G4_T1t200 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_kethoxal.T1t200.out'
	in_vitro_mRNA_T1t200 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_kethoxal.T1t200.out'
	in_vivo_G4_T0t20_out = read_icshape_out(in_vivo_G4_T0t20)
	in_vivo_mRNA_T0t20_out = read_icshape_out(in_vivo_mRNA_T0t20)
	in_vitro_G4_T1t200_out = read_icshape_out(in_vitro_G4_T1t200)
	in_vitro_mRNA_T1t200_out = read_icshape_out(in_vitro_mRNA_T1t200)

	fa_dict = read_fa('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/hg38/transcriptome/hg38_transcriptome.fa')

	save_dir = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/19-07-05/morestructure_overlap'

	for tx,start,end,mRNA_Gini1,G4_Gini1,mRNA_Gini2,G4_Gini2 in zip(df['tx'], df['start'], df['end'], df['mRNA gini_x'], df['G4 gini_x'],df['mRNA gini_y'], df['G4 gini_y']):
		seq = fa_dict[tx][start:end+extend]
		reactivity_ls3 = map(float, [-0.1 if i == 'NULL' else i for i in in_vivo_G4_T0t20_out[tx]['reactivity_ls'][start:end+extend]])
		reactivity_ls4 = map(float, [-0.1 if i == 'NULL' else i for i in in_vivo_mRNA_T0t20_out[tx]['reactivity_ls'][start:end+extend]])
		reactivity_ls1 = map(float, [-0.1 if i == 'NULL' else i for i in in_vitro_G4_T1t200_out[tx]['reactivity_ls'][start:end+extend]])
		reactivity_ls2 = map(float, [-0.1 if i == 'NULL' else i for i in in_vitro_mRNA_T1t200_out[tx]['reactivity_ls'][start:end+extend]])
		savefn = save_dir+'/'+'%s.%s-%s.pdf'%(tx, start, end)
		title = '%s:%s-%s, %.3f vs %.3f, %.3f vs %.3f'%(tx, start, end, mRNA_Gini1,G4_Gini1,mRNA_Gini2,G4_Gini2)
		ylim1 = max(reactivity_ls1+reactivity_ls2+reactivity_ls3+reactivity_ls4)
		# ylim2 = max(rt_ls1+rt_ls2)
		if tx != 'ENST00000356674': continue
		if start == 920: 
			print reactivity_ls1,reactivity_ls2,reactivity_ls3,reactivity_ls4
		else:
			continue
		plot_multiple_shape_track(seq=seq, reactivity_ls_ls=[reactivity_ls1, reactivity_ls2, reactivity_ls3, reactivity_ls4], highlight_base_ls=[], title=title, savefn=savefn, ylim_ls=[ylim1, ylim1, ylim1, ylim1], base='ATCG', set_non_G_0=1)

def plot_vivo_vitro_rt():
	# 19-04-03: in vivo reseq
	# vivo_G4_kethoxal_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4_kethoxal.rt'
	# vivo_mRNA_kethoxal_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_mRNA_kethoxal.rt'
	# vivo_G4_control_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4_control.rt'
	# vivo_mRNA_control_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_mRNA_control.rt'

	# use 1st batch
	vivo_G4_kethoxal_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4_kethoxal.rt'
	vivo_mRNA_kethoxal_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.rt'
	vivo_G4_control_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4_control.rt'
	vivo_mRNA_control_rt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_control.rt'

	vivo_G4_kethoxal_rt_dict = normRT_corr.read_rt(rt=vivo_G4_kethoxal_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)
	vivo_mRNA_kethoxal_rt_dict = normRT_corr.read_rt(rt=vivo_mRNA_kethoxal_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)
	vivo_G4_control_rt_dict = normRT_corr.read_rt(rt=vivo_G4_control_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)
	vivo_mRNA_control_rt_dict = normRT_corr.read_rt(rt=vivo_mRNA_control_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)

	vitro_G4_kethoxal_rt = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4_kethoxal.rt'
	vitro_mRNA_kethoxal_rt = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_mRNA_kethoxal.rt'
	vitro_G4_control_rt = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_control.rt'
	vitro_mRNA_control_rt = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_control.rt'
	vitro_G4_kethoxal_rt_dict = normRT_corr.read_rt(rt=vitro_G4_kethoxal_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)
	vitro_mRNA_kethoxal_rt_dict = normRT_corr.read_rt(rt=vitro_mRNA_kethoxal_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)
	vitro_G4_control_rt_dict = normRT_corr.read_rt(rt=vitro_G4_control_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)
	vitro_mRNA_control_rt_dict = normRT_corr.read_rt(rt=vitro_mRNA_control_rt, mode='raw', baseDensity_cutoff=0, gene_ls=None, RT_avg_threshold=0, split_id=1)

	fa_dict = read_fa('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/hg38/transcriptome/hg38_transcriptome.fa')
	extend = 20

	txt = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_kethoxal.T1t200.inRG4.txt'
	df = pd.read_csv(txt, header=None, sep='\t')
	df.columns = ['tx', 'start', 'end', 'gini', 'reactivity', 'seq']

	RT1_ls, RT2_ls = [], []
	RT1_gini_ls, RT2_gini_ls = [], []

	for tx,s,e in zip(df['tx'], map(int, df['start']), map(int, df['end'])):
		print tx,s,e
		start = int(s)
		end = int(e)

		vivo_G4_kethoxal_BD = map(float, vivo_G4_kethoxal_rt_dict[tx]['baseDensity'][start:end+extend])
		vivo_G4_control_BD = map(float, vivo_G4_control_rt_dict[tx]['baseDensity'][start:end+extend])
		vivo_mRNA_kethoxal_BD = map(float, vivo_mRNA_kethoxal_rt_dict[tx]['baseDensity'][start:end+extend])
		vivo_mRNA_control_BD = map(float, vivo_mRNA_control_rt_dict[tx]['baseDensity'][start:end+extend])

		vitro_G4_kethoxal_BD = map(float, vitro_G4_kethoxal_rt_dict[tx]['baseDensity'][start:end+extend])
		vitro_G4_control_BD = map(float, vitro_G4_control_rt_dict[tx]['baseDensity'][start:end+extend])
		vitro_mRNA_kethoxal_BD = map(float, vitro_mRNA_kethoxal_rt_dict[tx]['baseDensity'][start:end+extend])
		vitro_mRNA_control_BD = map(float, vitro_mRNA_control_rt_dict[tx]['baseDensity'][start:end+extend])

		seq = fa_dict[tx][start:end+extend]
		reactivity_ls_ls = [vivo_G4_kethoxal_BD, vivo_G4_control_BD, vivo_mRNA_kethoxal_BD, vivo_mRNA_control_BD,
							vitro_G4_kethoxal_BD, vitro_G4_control_BD, vitro_mRNA_kethoxal_BD, vitro_mRNA_control_BD]
		# ylim1 = max(reactivity_ls_ls[0]+reactivity_ls_ls[1]+reactivity_ls_ls[2]+reactivity_ls_ls[3])
		# ylim2 = max(reactivity_ls_ls[4]+reactivity_ls_ls[5]+reactivity_ls_ls[6]+reactivity_ls_ls[7])
		# ylim1 = max(reactivity_ls_ls[0]+reactivity_ls_ls[1])
		# ylim2 = max(reactivity_ls_ls[2]+reactivity_ls_ls[3])
		# ylim3 = max(reactivity_ls_ls[4]+reactivity_ls_ls[5])
		# ylim4 = max(reactivity_ls_ls[6]+reactivity_ls_ls[7])
		# ylim1 = max(reactivity_ls_ls[0])
		# ylim2 = max(reactivity_ls_ls[2])
		# ylim3 = max(reactivity_ls_ls[4])
		# ylim4 = max(reactivity_ls_ls[6])
		# ylim1 = max(reactivity_ls_ls[0])
		# ylim2 = max(reactivity_ls_ls[2])
		# ylim3 = max(reactivity_ls_ls[4])
		# ylim4 = max(reactivity_ls_ls[6])
		# RT1_ls.append(np.mean(RT1))
		# RT2_ls.append(np.mean(RT2))
		# RT1_gini_ls.append(gj.gini(RT1, mode='gini'))
		# RT2_gini_ls.append(gj.gini(RT2, mode='gini'))

		print tx,start,end,seq, reactivity_ls_ls
		# savefn = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-09-09_10_library_G4-glucose_time/PDS_rG4/%s:%s-%s.pdf'%(tx, start, end)
		# savefn = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/19-04-03/PDS_rG4/%s:%s-%s.png'%(tx, start, end)
		savefn = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/19-07-05/PDS_rG4/%s:%s-%s.pdf'%(tx, start, end)
		plot_multiple_shape_track(seq=seq, reactivity_ls_ls=reactivity_ls_ls, highlight_base_ls=[100],title='%s:%s-%s'%(tx, start, end), savefn=savefn, ylim_ls=None, set_non_G_0=0)

	# df['mean(RT1)'] = RT1_ls
	# df['mean(RT2)'] = RT2_ls
	# df['gini(RT1)'] = RT1_gini_ls
	# df['gini(RT2)'] = RT2_gini_ls
	# df.to_csv(txt.replace('.txt', '.rt.txt'), header=True, index=False, sep='\t')

def main():
	icshape_in_RG4('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/HeLa_kethoxal_vs_no-treat.shape.out', keep_G_only=0)
	icshape_in_RG4('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/HeLa_kethoxal_vs_no-treat.PDS.shape.out', keep_G_only=0)
	RG4_1='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/HeLa_kethoxal_vs_no-treat.shape.inRG4.txt'
	RG4_2='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/HeLa_kethoxal_vs_no-treat.PDS.shape.inRG4.txt'
	savefn_prefix = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/in_vivo_mRNA_rG4'
	RG4_compare(RG4_1=RG4_1, RG4_2=RG4_2, savefn_prefix=savefn_prefix)

	icshape_in_RG4('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/HeLa_kethoxal_vs_no-treat.invitro.shape.out', keep_G_only=0)
	icshape_in_RG4('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/HeLa_kethoxal_vs_no-treat.invitro.PDS.shape.out', keep_G_only=0)
	RG4_1='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/HeLa_kethoxal_vs_no-treat.invitro.shape.inRG4.txt'
	RG4_2='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/HeLa_kethoxal_vs_no-treat.invitro.PDS.shape.inRG4.txt'
	savefn_prefix = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/rG4/in_vitro_mRNA_rG4'
	RG4_compare(RG4_1=RG4_1, RG4_2=RG4_2, savefn_prefix=savefn_prefix)

if __name__ == '__main__':
	main()