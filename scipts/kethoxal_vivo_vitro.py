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
from scipy import stats

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
			rpkm = float(arr[2])
			reactivity_ls = arr[3:]
			out_dict[tx_id]['tx_id'] = tx_id
			out_dict[tx_id]['rpkm'] = rpkm
			out_dict[tx_id]['reactivity_ls'] = reactivity_ls
	gj.printFuncRun('read_icshape_out')
	return out_dict

def compare_structure_gini(out1=None, out2=None, condition1=None, condition2=None, save_dir=None, T=2, t=200):
	out1 = out1 if out1 is not None else '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/in_vivo_mRNA_kethoxal.T%st%s.out'%(T, t)
	out2 = out2 if out2 is not None else '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/in_vitro_mRNA_kethoxal.T%st%s.out'%(T, t)
	condition1 = condition1 if condition1 is not None else out1.split('/')[-1].split('.')[0].replace('_kethoxal', '')
	condition2 = condition2 if condition2 is not None else out2.split('/')[-1].split('.')[0].replace('_kethoxal', '')
	save_dir = save_dir if save_dir is not None else '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-08-08_16_library_invivo_invitro/gini'
	gj.printFuncRun('compare_structure_gini')
	gj.printFuncArgs()

	out_dict1 = read_icshape_out(out1)
	out_dict2 = read_icshape_out(out2)

	overlap_savefn = save_dir + '/' + '%s_%s.%s.txOverlap.png'%(condition1, condition2, out1.split('/')[-1].split('.')[-2])
	gj.venn3plot(mode='string',subsets_ls=[set(out_dict1.keys()), set(out_dict2.keys())],labels_ls=[condition1, condition2],title_str=None,save_fn=overlap_savefn,axis=None)

	overlap_tx = set(out_dict1.keys()) & set(out_dict2.keys())
	null_pct_ls = [0.2, 0.4, 0.6, 0.8, 0.9, 1]
	gini_savefn = save_dir + '/' + '%s-%s.txt'%(out1.split('/')[-1], out2.split('/')[-1])
	SAVEFN = open(gini_savefn, 'w')
	header_ls = ['tx', 'null_pct_cutoff', 'vivo', 'vitro', 'null_pct1', 'null_pct2']
	print >>SAVEFN, '\t'.join(header_ls)
	for null_pct in null_pct_ls:
		out1_gini_ls = [gj.gini(out_dict1[tx]['reactivity_ls'], mode='gini',null_pct=null_pct) for tx in overlap_tx]
		out2_gini_ls = [gj.gini(out_dict2[tx]['reactivity_ls'], mode='gini',null_pct=null_pct) for tx in overlap_tx]

		out1_gini_ls_filter = []
		out2_gini_ls_filter = []
		overlap_tx_filter = []
		for i,j,tx in zip(out1_gini_ls, out2_gini_ls, overlap_tx):
			if float(i) >= 0 and float(j) >= 0:
				out1_gini_ls_filter.append(i)
				out2_gini_ls_filter.append(j)
				overlap_tx_filter.append(tx)

				out1_tx_null_pct = len(['NULL' for n in out_dict1[tx]['reactivity_ls'] if n == 'NULL']) / float(len(out_dict1[tx]['reactivity_ls']))
				out2_tx_null_pct = len(['NULL' for n in out_dict2[tx]['reactivity_ls'] if n == 'NULL']) / float(len(out_dict2[tx]['reactivity_ls']))

				print >> SAVEFN, '\t'.join(map(str, [tx, null_pct, i, j, out1_tx_null_pct, out2_tx_null_pct]))
	SAVEFN.close()

	gj.printFuncRun('compare_structure_gini')

def plot_gini_compare(save_dir=None):
	save_dir = save_dir if save_dir is not None else '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/result/16-08-08_16_library_invivo_invitro/gini'
	fn_ls = os.listdir(save_dir)
	T_t_ls = ['T%st%s'%(i,j) for i in [0,1,2] for j in [0,20,200]]
	fn_ls = [i for i in fn_ls if i.endswith('.txt') and i.split('.')[1] in T_t_ls]
	print fn_ls

	df_ls = []
	for fn in fn_ls:
		fn_path = save_dir + '/' + fn
		print "process: %s"%(fn_path)
		df = pd.read_csv(fn_path, header=0, sep='\t')
		condition = fn.split('.')[0].split('_')[2]
		T_t_cutoff = fn.split('.')[1]
		df['condition'] = condition
		df['T_t_cutoff'] = T_t_cutoff

		print df.head()

		df_ls.append(df)
	df_all = pd.concat(df_ls)
	savefn = save_dir + '/' + 'vivo_vitro_gini_all.txt'
	df_all.to_csv(savefn, index=False, header=True, sep='\t')
	print df_all.head()

	df_all_melt = pd.melt(df_all, id_vars=['tx', 'null_pct_cutoff', 'condition', 'T_t_cutoff'], value_vars=['vivo', 'vitro'], var_name='vivo/vitro', value_name='Gini')
	df_all_melt['condition'] = ['%s,%s'%(i,j) for i,j in zip(df_all_melt['condition'], df_all_melt['vivo/vitro'])]
	print df_all_melt.head()
	df_all_melt.sort_values(by=['condition', 'T_t_cutoff'], inplace=True)
	df_all_melt = df_all_melt[df_all_melt['null_pct_cutoff'].isin([0.4, 0.6])]
	df_all_melt = df_all_melt[df_all_melt['T_t_cutoff'].isin(['T0t0', 'T0t20', 'T1t20'])]
	df_all_melt.to_csv(savefn.replace('.txt', '.plot.txt'), index=False, header=True, sep='\t')
	g = sns.FacetGrid(data=df_all_melt, row='T_t_cutoff', col='null_pct_cutoff', sharey=True, margin_titles=True)
	g = g.map(sns.boxplot, 'condition', 'Gini', )
	g.set_xticklabels(rotation=90)
	g.set(ylim=(0.5,1))
	g.savefig(savefn.replace('.txt', '.png'))
	plt.close()

	g = sns.FacetGrid(data=df_all_melt, row='T_t_cutoff', col='null_pct_cutoff', sharey=True, margin_titles=True)
	g = g.map(sns.countplot, 'condition',)
	g.set_xticklabels(rotation=90)
	g.savefig(savefn.replace('.txt', '.count.png'))
	plt.close()

	df_select = df_all_melt[(df_all_melt['null_pct_cutoff']==0.4) & (df_all_melt['T_t_cutoff']=='T1t20') & df_all_melt['condition'].isin(['mRNA,vivo', 'mRNA,vitro'])]
	print df_select

	fig,ax=plt.subplots(figsize=(4,6))
	sns.boxplot(x='condition', y='Gini', data=df_select, ax=ax)
	plt.tight_layout()
	ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45)
	plt.savefig(savefn.replace('.txt', '.select.png'))
	plt.close()
	
	condition_ls = []
	condition_gini_ls = []
	for i in df_select['condition'].value_counts().keys():
		print i
		df_select[df_select['condition']==i]
		condition_ls.append(i)
		condition_gini_ls.append(list(df_select[df_select['condition']==i]['Gini']))
		savefn = savefn.replace('.txt', '.select2.png')
	gj.cumulate_dist_plot(ls_ls=condition_gini_ls,ls_ls_label=condition_ls,bins=40,title=None,ax=None,savefn=savefn,xlabel=None,ylabel=None,add_vline=None,add_hline=None,log2transform=0)
	
	stat,pval = stats.ks_2samp(condition_gini_ls[0],condition_gini_ls[1])
	print pval

def main():
	# compare_structure_gini()

	# """
	T_ls = [1, ]
	t_ls = [20, ]
	save_dir = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/vivo_vitro_gini'
	for T,t in zip(T_ls, t_ls):
		compare_structure_gini(T=T, t=t, save_dir=save_dir)
	# """

	# plot_gini_compare()

if __name__ == '__main__':
	main()
