import gj
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
plt.rcParams["font.family"] = "Helvetica"
import sys, os
from nested_dict import nested_dict
import pandas as pd
import numpy as np
from pyfasta import Fasta
from scipy import stats

def read_fa(fa='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'):
	gj.printFuncRun('read_fa')
	gj.printFuncArgs()
	fa_dict1 = Fasta(fa, key_fn=lambda key:key.split("\t")[0])
	fa_dict = {i.split()[0].split('.')[0]:j[0:] for i,j in fa_dict1.items()}
	print fa_dict.keys()[0:3]
	gj.printFuncRun('read_fa')
	return fa_dict

def read_bed(bed=None):
	if bed is None:
		bed = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10.transCoor.bed'
	bed_dict = {}
	with open(bed, 'r') as BED:
		for line in BED:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			arr = line.split('\t')
			tx_id = arr[0].split('.')[0]
			gene_name = arr[1].split('=')[0]
			bed_dict[tx_id] = gene_name
	print "bed",bed,bed_dict.keys()[0:3]
	return bed_dict

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
			out_dict[tx_id]['rpkm'] = rpkm
			out_dict[tx_id]['reactivity_ls'] = reactivity_ls
	gj.printFuncRun('read_icshape_out')
	return out_dict

def compare_corr(out1=None, out2=None, savefn=None, label=None):
	fa_dict = read_fa('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa')

	label1, label2 = label.split(':')
	out1 = out1 if out1 is not None else '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.out'
	out2 = out2 if out2 is not None else '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/data/icSHAPE/icshape.out'

	out_dict1 = read_icshape_out(out1)
	out_dict2 = read_icshape_out(out2)

	labels_ls = [out1.split('/')[-1], out2.split('/')[-1]]
	# gj.venn3plot(mode='string',subsets_ls=[set(out_dict1.keys()), set(out_dict2.keys())],labels_ls=labels_ls,title_str=None,save_fn=savefn.replace('.txt', '.tx_overlap.png'),axis=None)

	tx_overlap_ls = set(out_dict1.keys()) & set(out_dict2.keys())

	SAVEFN = open(savefn, 'w')
	header_ls = ['tx', 'pos', 'base',label1, label2]
	print >>SAVEFN, '\t'.join(header_ls)
	for tx in tx_overlap_ls:
		if len(out_dict1[tx]['reactivity_ls']) != len(out_dict2[tx]['reactivity_ls']):
			print "not equal len: %s, %s, %s"%(tx, len(out_dict1[tx]['reactivity_ls']), len(out_dict2[tx]['reactivity_ls']))
			continue
		for n,(i,j) in enumerate(zip(out_dict1[tx]['reactivity_ls'], out_dict2[tx]['reactivity_ls'])):
			if i != 'NULL' and j != 'NULL':
				print >>SAVEFN, '\t'.join(map(str, [tx, n, fa_dict[tx][n:n+1], i, j]))
	SAVEFN.close()

	df = pd.read_csv(savefn, header=0, sep='\t')
	# print df.head()
	# fig,ax = plt.subplots()
	# sns.lmplot(x=label1, y=label2, data=df, hue='base', size=4, col='base', col_wrap=2, scatter_kws={"s": 5})
	# #df.plot(kind='scatter', x='reactivity1', y='reactivity2', ax=ax)
	# plt.savefig(savefn.replace('.txt', '.png'))
	# plt.close()

	"""
	for b in ['A', 'T', 'C', 'G']:
		df_base = df[df['base']==b]
		print df_base.head()
		g = sns.jointplot(x=label1, y=label2, data=df_base, kind="kde", scatter_kws={"s": 5}, xlim=(-0.05,1.05), ylim=(-0.05,1.05), color='#8172B2')
		g.ax_joint.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
		sns.regplot(x=label1, y=label2, data=df_base, scatter=False, ax=g.ax_joint, color='#8172B2')
		plt.savefig(savefn.replace('.txt', '.%s.png'%(b)))
		plt.close()
	"""

	# tx_corr_savefn = savefn.replace('.txt', '.list.txt')
	# TX_CORR_SAVEFN = open(tx_corr_savefn, 'w')
	bed_dict = read_bed('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10.transCoor.bed')

	tx_ls = df['tx'].value_counts().keys()
	tx_corr_ls = []
	tx_base_corr_ls_ls = [[], [], [], []]
	for tx in tx_ls:
		df_tx = df[df['tx']==tx]
		if df_tx.shape[0] <= 10:
			continue
		else:
			df_tx_corr = df_tx[[label1, label2]].corr()
			tx_corr_ls.append(df_tx_corr.iloc[0,1])
		
		for n,base in enumerate(['A', 'T', 'C', 'G']):
			df_tx_base = df_tx[df_tx['base']==base]
			if df_tx_base.shape[0] <= 10:
				pass
			else:
				df_tx_base_corr = df_tx_base[[label1, label2]].corr()
				if np.isnan(df_tx_base_corr.iloc[0,1]):
					pass
				else:
					tx_base_corr_ls_ls[n].append(df_tx_base_corr.iloc[0,1])
					base_pct = fa_dict[tx][0:].count(base) / float(len(fa_dict[tx][0:]))
					GC_content = (fa_dict[tx][0:].count('C') + fa_dict[tx][0:].count('G')) / float(len(fa_dict[tx][0:]))
					# print >>TX_CORR_SAVEFN, '\t'.join(map(str, [tx, bed_dict[tx], base, base_pct, df_tx_base.shape[0], df_tx_base_corr.iloc[0,1], GC_content]))
	# TX_CORR_SAVEFN.close()

	return [tx_base_corr_ls_ls[-1]]+[tx_corr_ls]

def tx_corr_scatter(txt=None, tx='ENSMUST00000201155'):
	if txt is None:
		txt = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/F2b.kethoxal-notreat.icshape.corr.txt'
	df = pd.read_csv(txt, header=0, sep='\t')
	print df.head()

	df = df[df['tx']==tx]
	df = df[df['base']=='G']

	fig,ax = plt.subplots(figsize=(6,6))
	df.plot(kind='scatter', x='kethoxalseq', y='icSHAPE', ax=ax)
	R,p = stats.pearsonr(df['kethoxalseq'], df['icSHAPE'])
	plt.title('R=%.3f, p=%s'%(R, p))
	ax.set_aspect('equal')
	ax.set_xlim(-0.05, 1.05)
	ax.set_ylim(-0.05, 1.05)
	plt.tight_layout()
	savefn = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/F2b.%s.corr.scatter.pdf'%(tx)
	plt.savefig(savefn)
	plt.close()


def main():
	notreat = compare_corr(out1='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.out', savefn='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/F2b.kethoxal-notreat.icshape.corr.txt', label='kethoxalseq:icSHAPE')
	gj.cumulate_dist_plot(ls_ls=[notreat[-2]],ls_ls_label=['kethoxal vs icshape', ],bins=40,title=None,ax=None,savefn='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/F2b.kethoxal_icshape_corr.pdf',xlabel=None,ylabel=None,add_vline=None,add_hline=None,log2transform=0)
	tx_corr_scatter()

if __name__ == '__main__':
	main()