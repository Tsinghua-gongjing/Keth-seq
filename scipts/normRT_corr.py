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
from scipy import stats

def read_rpkm_txt(txt, min_val=-1):
	gj.printFuncRun('read_rpkm_txt')
	gj.printFuncArgs()
	val_dict = nested_dict()
	gene_ls = []
	with open(txt, 'r') as TXT:
		for line in TXT:
			line = line.strip()
			if not line or line.startswith('#'): continue
			arr = line.split('\t')
			val_dict[arr[0]] = float(arr[4])
			gene_ls.append(arr[0])
	gj.printFuncRun('read_rpkm_txt')
	return val_dict,gene_ls

def read_fa(fa='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'):
	gj.printFuncRun('read_fa')
	gj.printFuncArgs()
	fa_dict = Fasta(fa, key_fn=lambda key:key.split("\t")[0])
	print fa_dict.keys()[0:3]
	gj.printFuncRun('read_fa')
	return fa_dict

# def read_rt(rt=None):
# 	gj.printFuncRun('read_rt')
# 	gj.printFuncArgs()
# 	rt_dict = nested_dict()
# 	with open(rt, 'r') as RT:
# 		for line in RT:
# 			line = line.strip()
# 			if not line or line.startswith('#'): continue
# 			arr = line.split('\t')
# 			tx_id = arr[0]
# 			rt_dict[tx_id]['tx_id'] = tx_id
# 			rt_dict[tx_id]['length'] = int(arr[1])
# 			rt_dict[tx_id][arr[2]] = arr[5:] 
# 	#gj.print_dict_item0(rt_dict)
# 	gj.printFuncRun('read_rt')
# 	return rt_dict.to_dict()

def read_rt(rt=None, mode='normalized', baseDensity_cutoff=2, gene_ls=None, RT_avg_threshold=2, split_id=0):
	print "read RT: %s"%(rt)
	rt_dict = nested_dict()
	
	with open(rt, 'r') as RT:
		n = 0
		for line in RT:
			line = line.strip()
			if not line or line.startswith('#'): continue
			n += 1
			arr = line.split('\t')
			tx_id = arr[0].split('.')[0] if split_id else arr[0]
			if mode == 'raw':
				rt_dict[tx_id]['tx_id'] = tx_id
				rt_dict[tx_id]['length'] = int(arr[1])
				rt_dict[tx_id]['rpkm'] = float(arr[2])
				if n%2 == 1:
					rt_dict[tx_id]['baseDensity'] = arr[4:]
				if n%2 == 0:
					rt_dict[tx_id]['RTstop'] = arr[4:]
			if mode == 'normalized':
				rt_dict[tx_id]['tx_id'] = tx_id
				rt_dict[tx_id]['length'] = int(arr[1])
				rt_dict[tx_id][arr[2]] = arr[5:]
				rt_dict[tx_id][arr[2]+'_scalefactor'] = arr[4]

	if gene_ls is not None:
		gene_ls_dict = {i:0 for i in gene_ls}
		for i,j in rt_dict.items():
			if not gene_ls_dict.has_key(i):
				rt_dict.pop(i)
	for i,j in rt_dict.items():
		if not j['RTstop']:
			rt_dict.pop(i)
	if mode == 'raw':
		for i,j in rt_dict.items():
			BD_avg = np.mean(map(float, j['baseDensity']))
			RT_avg = np.mean(map(float, j['RTstop']))
			if BD_avg < baseDensity_cutoff or RT_avg < RT_avg_threshold:
				rt_dict.pop(i)
	return rt_dict.to_dict()

def compare_RT(rt1=None, rt2=None, savefn=None):
	if rt1 is None:
		# rt1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1C_S1_L006_R1_001.normalized.rt'
		# rt1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.rt'
		rt1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.rt'
		# rt1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1C_S1_L006_R1_001.rt'
	if rt2 is None:
		# rt2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2C_S2_L006_R1_001.normalized.rt'
		# rt2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.rt'
		rt2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON2_S6_L006_R1_001.rt'
		# rt2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2C_S2_L006_R1_001.rt'
	if savefn is None:
		# savefn = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-11-14_7_library_total_Kethoxal_remove/kethoxal_rep_raw_RT_corr.T10t1.txt'
		savefn = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-11-14_7_library_total_Kethoxal_remove/notreat_rep_raw_RT_corr.T1t1.txt'
	rt_dict1 = read_rt(rt1, mode='raw', baseDensity_cutoff=1, RT_avg_threshold=1)
	rt_dict2 = read_rt(rt2, mode='raw', baseDensity_cutoff=1, RT_avg_threshold=1)
	fa_dict = read_fa()

	rt1_tx_ls = [i for i in rt_dict1 if rt_dict1[i].has_key('RTstop')]
	rt2_tx_ls = [i for i in rt_dict2 if rt_dict2[i].has_key('RTstop')]

	rt1_rt2_common_tx_ls = set(rt1_tx_ls) & set(rt2_tx_ls)

	SAVEFN = open(savefn, 'w')
	head_ls = ['tx', 'tx_index', 'base', 'RTstop1', 'RTstop2']
	print >>SAVEFN, '\t'.join(head_ls) 
	for tx in rt1_rt2_common_tx_ls:
		if len(rt_dict1[tx]['RTstop']) != len(rt_dict2[tx]['RTstop']):
			continue
		if len(rt_dict1[tx]['RTstop']) != len(fa_dict[tx][0:]):
			continue
		for n,(RTstop1,RTstop2,base) in enumerate(zip(rt_dict1[tx]['RTstop'], rt_dict2[tx]['RTstop'], fa_dict[tx][0:])):
			print >>SAVEFN, '\t'.join(map(str, [tx, n, base, RTstop1, RTstop2]))
	SAVEFN.close()

	print "rt1: %s"%(rt1)
	print "   - tx num: %s, has RT: %s"%(len(rt_dict1), len(rt1_tx_ls))
	print "rt2: %s"%(rt2)
	print "   - tx num: %s, has RT: %s"%(len(rt_dict2), len(rt2_tx_ls))
	print "\nCommon tx num: %s"%(len(rt1_rt2_common_tx_ls))

def compare_RT_plot(compare_RT_txt=None):
	if compare_RT_txt is None:
		compare_RT_txt = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-11-14_7_library_total_Kethoxal_remove/notreat_rep_RT_corr.txt'
	RT_rep_rpkm_dict = {'notreat_rep_RT_corr.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.rpkm',
												   '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON2_S6_L006_R1_001.rpkm'],
						'kethoxal_rep_RT_corr.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.rpkm',
													'/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.rpkm'],
						'control_rep_RT_corr.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1C_S1_L006_R1_001.rpkm',
												   '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2C_S2_L006_R1_001.rpkm'],
						'notreat_rep_raw_RT_corr.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.rpkm',
												   '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON2_S6_L006_R1_001.rpkm'],
						'notreat_rep_raw_RT_corr.T0t0.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.rpkm',
												   '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON2_S6_L006_R1_001.rpkm'],
						'notreat_rep_raw_RT_corr.T0t1.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.rpkm',
												   '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON2_S6_L006_R1_001.rpkm'],
						'notreat_rep_raw_RT_corr.T1t1.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.rpkm',
												   '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON2_S6_L006_R1_001.rpkm'],
						'notreat_rep_raw_RT_corr.T10t1.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.rpkm',
												   '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON2_S6_L006_R1_001.rpkm'],
						'kethoxal_rep_raw_RT_corr.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.rpkm',
													'/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.rpkm'],
						'kethoxal_rep_raw_RT_corr.T0t0.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.rpkm',
													'/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.rpkm'],
						'kethoxal_rep_raw_RT_corr.T0t1.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.rpkm',
													'/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.rpkm'],
						'kethoxal_rep_raw_RT_corr.T1t1.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.rpkm',
													'/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.rpkm'],
						'kethoxal_rep_raw_RT_corr.T10t1.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.rpkm',
													'/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.rpkm'],
						'control_rep_raw_RT_corr.txt':['/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1C_S1_L006_R1_001.rpkm',
												   '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2C_S2_L006_R1_001.rpkm']}
	rpkm_ls = RT_rep_rpkm_dict[compare_RT_txt.split('/')[-1]]
	rpkm_dict_ls = []
	for rpkm in rpkm_ls:
		val_dict,gene_ls = read_rpkm_txt(rpkm)
		rpkm_dict_ls.append(val_dict)

	df = pd.read_csv(compare_RT_txt, header=0, sep='\t')
	df['allGt0'] = [0 if i==0 and j==0 else 1 for i,j in zip(list(df['RTstop1']), list(df['RTstop2']))]
	# df['allGt0'] = [1 if i > 0 and j > 0 else 0 for i,j in zip(list(df['RTstop1']), list(df['RTstop2']))]
	df = df[df['allGt0'] == 1]
	df_rpkm = pd.DataFrame.from_dict({'rep1':[np.log2(rpkm_dict_ls[0][i]) for i in set(list(df['tx']))], 'rep2':[np.log2(rpkm_dict_ls[1][i]) for i in set(list(df['tx']))]})
	fig, ax = plt.subplots(figsize=(3,5))
	df_rpkm.plot(kind='box', ax=ax)
	plt.tight_layout()
	plt.savefig(compare_RT_txt.replace('.txt', '.rpkm.png'))
	plt.close()
	rpkm_min_val = 0
	df['rpkm'] = [1 if rpkm_dict_ls[0][i] > rpkm_min_val and rpkm_dict_ls[0][i] > rpkm_min_val else 0 for i in list(df['tx'])]
	df['rpkm1'] = [rpkm_dict_ls[0][i] for i in list(df['tx'])]
	df['rpkm2'] = [rpkm_dict_ls[1][i] for i in list(df['tx'])]
	df = df[df['rpkm']==1]
	df.replace(0, 0.25, inplace=True)
	df['log2(RTstop1)'] = np.log2(df['RTstop1'])
	df['log2(RTstop2)'] = np.log2(df['RTstop2'])
	print df.head()
	# df = df[df['base']!='G']
	print df['base'].value_counts()
	print "tx num: %s, uniq: %s"%(len(df['tx']), len(set(df['tx'])))
	print df.describe()

	r,p = stats.pearsonr(x=df['RTstop1'],y=df['RTstop2']) # return (r, p)
	print p,r 

	fig,ax = plt.subplots(figsize=(6,6))
	# base_ls = df['base'].value_counts().keys()
	base_ls = ['A', 'T', 'C', 'G']
	color_ls = ['#4C72B0', '#55A868', '#C44E52', '#8172B2']
	# color_ls = ['#d02e30', '#57a455', '#e9e847', 'grey']
	color_ls_dict = dict(zip(base_ls, color_ls))
	df['color'] = [color_ls_dict[i] for i in df['base']]
	"""
	for n,base in enumerate(base_ls):
		df_base = df[df['base']==base]
		df_base.plot(kind='scatter', x='log2(RTstop1)', y='log2(RTstop2)', color=color_ls[n], ax=ax, label=base, alpha=0.5, s=0.8)
	"""
	"""
	for x,y,b in zip(df['log2(RTstop1)'], df['log2(RTstop2)'], df['base']):
		ax.scatter(x,y,c=color_ls_dict[b],alpha=0.3, s=0.8)
	"""

	#sns.pointplot(x='log2(RTstop1)', y='log2(RTstop2)', data=df.iloc[0:10000,:], hue='base', linestyles='', ax=ax, ci=None)
	ax.scatter(df['log2(RTstop1)'], df['log2(RTstop2)'], c=df['color'], alpha=0.5, s=0.8)
	# g = sns.FacetGrid(df, col='base', hue='base', sharex=True, sharey=True)
	# g.map(sns.kdeplot, 'log2(RTstop1)', 'log2(RTstop2)')

	savefn = compare_RT_txt.replace('.txt', '.png')
	plt.title('R2=%s, pvalue=%s'%(r, p))
	ax.axis('square')
	ax.set_xlim(-3.5, 15.5)
	ax.set_ylim(-3.5, 15.5)
	# ax.legend_.remove()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	plt.xticks([-3, 0, 3, 6, 9, 12, 15])
	plt.yticks([-3, 0, 3, 6, 9, 12, 15])
	plt.tight_layout()
	plt.savefig(savefn, dpi=220)
	plt.close()

def main():
	# compare_RT()
	# compare_RT_plot('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-11-14_7_library_total_Kethoxal_remove/kethoxal_rep_RT_corr.txt')
	# compare_RT_plot('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-11-14_7_library_total_Kethoxal_remove/control_rep_RT_corr.txt')
	# compare_RT_plot('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-11-14_7_library_total_Kethoxal_remove/notreat_rep_raw_RT_corr.T0t0.txt')
	compare_RT_plot('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/F2a.kethoxal_rep_RT_corr.txt')

if __name__ == '__main__':
	main()