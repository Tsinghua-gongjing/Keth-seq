import gj
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
import sys
from nested_dict import nested_dict
import pandas as pd
import numpy as np
from pyfasta import Fasta
import subprocess

def read_fa(fa='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'):
	gj.printFuncRun('read_fa')
	gj.printFuncArgs()
	fa_dict = Fasta(fa, key_fn=lambda key:key.split("\t")[0])
	print fa_dict.keys()[0:3]
	gj.printFuncRun('read_fa')
	return fa_dict

def base_complementary(base='A'):
	return {'A':'T','T':'A','C':'G','G':'C'}[base]

def read_bed(bed=None, fa=None):
	gj.printFuncRun('read_bed')
	gj.printFuncArgs()
	base_dict = nested_dict(1, int)
	fa_dict = read_fa(fa)
	with open(bed, 'r') as BED:
		for n,line in enumerate(BED):
			line = line.strip()
			if not line or line.startswith('#'): continue
			if n%1000000 == 0: print "process: %s"%(n)
			arr = line.split('\t')
			tx_id = arr[0]
			tx_start = int(arr[1])
			tx_end = int(arr[2])
			strand = arr[5]
			if strand == "+":
				base = fa_dict[tx_id][tx_start-1]
			elif strand == "-":
				# base = fa_dict[tx_id][tx_end]
				# base = base_complementary(base)
				continue
			else:
				print "unknown strand: %s"%(strand)
				sys.exit()
			base_dict[base] += 1
	print base_dict

	bed_base_txt = bed.replace('bed','base.txt')
	with open(bed_base_txt,'w') as TXT:
		for i,j in base_dict.items():
			print >>TXT,i+'\t'+str(j)

	gj.printFuncRun('read_bed')

def base_ratio_plot(dir='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-08-14', base_txt_ls=None, condition_ls=None, savefn=None):
	gj.printFuncRun('base_ratio_plot')
	# gj.printFuncArgs()

	# base_txt_ls = ['CHe-XC-CON1_S5_L006_R1_001.map.plus.base.txt','CHe-XC-CON2_S6_L006_R1_001.map.plus.base.txt', 'CHe-XC-M1C_S1_L006_R1_001.map.plus.base.txt','CHe-XC-M2C_S2_L006_R1_001.map.plus.base.txt','CHe-XC-M1K_S3_L006_R1_001.map.plus.base.txt',
	# 				'CHe-XC-M2K_S4_L006_R1_001.map.plus.base.txt',]
	# condition_ls = ['No treat1', 'No treat2', 'Control rep1','Control rep2','Kethoxal rep1','Kethoxal rep2']

	# base_txt_ls = ['C1_L4_A021.R1.map.plus.base.txt', 'C2_L4_A022.R1.map.plus.base.txt', 'K1_L4_A024.R1.map.plus.base.txt', 'K2_L4_A025.R1.map.plus.base.txt']
	# condition_ls = ['mRNA_control', 'G4_control', 'mRNA_kethoxal', 'G4_kethoxal']

	# base_txt_ls = ['NpdNP_combined_R2.map.plus.base.txt', 'NpdP_combined_R2.map.plus.base.txt', 'pdNP_combined_R2.map.plus.base.txt', 'pdP_combined_R2.map.plus.base.txt']
	# condition_ls = ['rG4 kethxoal (-pulldown, +PDS)', 'kethxoal (-pulldown, -PDS)', 'kethoxal (+pulldown, -PDS)', 'rG4 kethoxal (+pulldown, +PDS)']

	# base_txt_ls = ['YSXK1_combined_R2.map.plus.base.txt', 'YSXK_combined_R2.map.plus.base.txt', 'YSXNK1_combined_R2.map.plus.base.txt', 'YSXNK_combined_R2.map.plus.base.txt']
	# condition_ls = ['YSXK1', 'YSXK', 'YSXNK1', 'YSXNK']

	# base_txt_ls = ['CHe-Weng-C-1_S1_L006_R1_001.map.plus.base.txt', 'CHe-Weng-C-2_S2_L006_R1_001.map.plus.base.txt', 'CHe-Weng-K-1_S4_L006_R1_001.map.plus.base.txt', 'CHe-Weng-K-2_S5_L006_R1_001.map.plus.base.txt']
	# condition_ls = ['mRNA_control', 'G4_control', 'mRNA_kethoxal', 'G4_kethoxal']

	# base_txt_ls = ['CHe-XC-CON1_S5_L006_R1_001.map.plus.base.txt', 'CHe-Weng-T-1_S7_L006_R1_001.map.plus.base.txt',
	# 				'CHe-Weng-T-2_S8_L006_R1_001.map.plus.base.txt', 'CHe-Weng-T-3_S9_L006_R1_001.map.plus.base.txt', 'CHe-Weng-T-4_S10_L006_R1_001.map.plus.base.txt']
	# condition_ls = ['Control', 'Label 1 min', 'Label 2.5 min', 'Label 5 min', 'Label 10 min']

	# base_txt_ls = ['G1_combined_R2.map.plus.base.txt', 'G2_combined_R2.map.plus.base.txt', 'K1_combined_R2.map.plus.base.txt',
	# 			   'K2_combined_R2.map.plus.base.txt', 'M1_combined_R2.map.plus.base.txt', 'M2_combined_R2.map.plus.base.txt',
	# 			   'mG1_combined_R2.map.plus.base.txt', 'mG2_combined_R2.map.plus.base.txt', 'mK1_combined_R2.map.plus.base.txt',
	# 			   'mK2_combined_R2.map.plus.base.txt', 'mM1_combined_R2.map.plus.base.txt', 'mM2_combined_R2.map.plus.base.txt']
	# condition_ls = ['Glyxal1(total)','Glyxal2(total)','Kethoxal1(total)', 'Kethoxal2(total)', 'Methylgyxal1(total)', 'Methylgyxal2(total)',
	# 				'Glyxal1(mRNA)','Glyxal2(mRNA)','Kethoxal1(mRNA)', 'Kethoxal2(mRNA)', 'Methylgyxal1(mRNA)', 'Methylgyxal2(mRNA)',]

	df_ls = []
	for base_txt,condition in zip(base_txt_ls,condition_ls):
		df = pd.read_csv(dir+'/'+base_txt, sep='\t', header=None)
		df.columns = ['base',condition]
		df.index = df['base']
		print df 
		df_ls.append(df)
		
	df_all = pd.concat(df_ls, axis=1)
	df_all = df_all[condition_ls].T
	df_all = df_all[['A', 'T', 'C', 'G']]
	print df_all
	fig,ax = plt.subplots(figsize=(12, 6))
	df_all.plot.barh(stacked=True, ax=ax)
	# ax.set_title('Base count at stop site(called from sam) across different samples')
	plt.tight_layout()
	plt.savefig(savefn)
	plt.close()

	print df_all.sum(axis=1)
	df_all.sum(axis=1).to_csv(savefn.replace('.pdf', '.sum.txt'))

	df_all_normalize = df_all.div(df_all.sum(axis=1), axis=0)
	print df_all_normalize

	fig,ax = plt.subplots(figsize=(12, 6))
	df_all_normalize.plot.barh(stacked=True,ax=ax).legend(loc='center left')
	# df_all.plot.barh(stacked=True,ax=ax).legend(loc='center left', bbox_to_anchor=(1, 0.5))
	# ax.set_title('Base ratio at stop site(called from sam) across different samples')
	ax.set_xlim(0,1)
	plt.tight_layout()
	plt.savefig(savefn.replace('.pdf', '.ratio.pdf'))
	plt.close()

	gj.printFuncRun('base_ratio_plot')

if __name__ == "__main__":
	""" python sam_RT_base_ratio.py .sam transcriptome.fa
	"""
	# sam = sys.argv[1]
	# fa = sys.argv[2]
	# bed = sam.replace('sam','bed')
	# subprocess.call(["sam2bed %s %s"%(sam, bed)],shell=True)
	# bed = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.map.plus.bed'
	# bed = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-T-1_S7_L006_R1_001.map.plus.bed'
	# fa = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'

	# bed = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.map.plus.bed'
	# fa = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'
	# read_bed(bed=bed, fa=fa)

	bed = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/test.map.plus.bed'
	fa = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'
	read_bed(bed=bed, fa=fa)

	# base_txt_ls = ['CHe-XC-CON1_S5_L006_R1_001.map.plus.base.txt', 'CHe-Weng-T-1_S7_L006_R1_001.map.plus.base.txt',
	# 				'CHe-Weng-T-2_S8_L006_R1_001.map.plus.base.txt', 'CHe-Weng-T-3_S9_L006_R1_001.map.plus.base.txt', 'CHe-Weng-T-4_S10_L006_R1_001.map.plus.base.txt']
	# condition_ls = ['Control', 'Label 1 min', 'Label 2.5 min', 'Label 5 min', 'Label 10 min']
	# base_ratio_plot(dir='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time', base_txt_ls=base_txt_ls, condition_ls=condition_ls, savefn='../results/SF5.time_points.RT.count.pdf')


	# base_txt_ls = ['CHe-XC-CON1_S5_L006_R1_001.map.plus.base.txt','CHe-XC-CON2_S6_L006_R1_001.map.plus.base.txt', 'CHe-XC-M1C_S1_L006_R1_001.map.plus.base.txt','CHe-XC-M2C_S2_L006_R1_001.map.plus.base.txt','CHe-XC-M1K_S3_L006_R1_001.map.plus.base.txt',
	# 				'CHe-XC-M2K_S4_L006_R1_001.map.plus.base.txt',]
	# condition_ls = ['No treat1', 'No treat2', 'Remove rep1','Remove rep2','Kethoxal rep1','Kethoxal rep2']
	# base_ratio_plot(dir='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove', base_txt_ls=base_txt_ls, condition_ls=condition_ls, savefn='../results/SF8b.rep.RT.count.pdf')
