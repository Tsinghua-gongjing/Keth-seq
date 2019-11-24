from nested_dict import nested_dict
import gj
import pandas as pd, numpy as np
from sklearn.metrics import roc_curve, auc
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")
sns.set_context("poster")
plt.rcParams["font.family"] = "Helvetica"
from pyfasta import Fasta
import random 

def read_dot(dot):
	seq_dict = nested_dict(2, list)
	with open(dot, 'r') as DOT:
		for line in DOT:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			if line.startswith('>'):
				seq_id = line.split('\t')[0].replace('>', '')
				seq = DOT.next().strip()
				dotstr = DOT.next().strip()
				seq_dict[seq_id]['seq_id'] = line
				seq_dict[seq_id]['seq'] = seq
				seq_dict[seq_id]['dotstr'] = dotstr
	# print seq_dict.keys()
	return seq_dict.to_dict()

def simulate_dot(dot, ss_min=0.8, ss_max=1, ds_min=0, ds_max=0.2, savefn=None):
	random.seed(1234)
	seq_dict = read_dot(dot)
	savefn = savefn.replace('.txt', '.ss_min%s.ss_max%s.ds_min%s.ds_max%s.txt'%(ss_min, ss_max, ds_min, ds_max))
	with open(savefn, 'w') as SAVEFN:
		for i,j in seq_dict.items():
			dot = j['dotstr']
			ss = j['dotstr'].count('.') / float(len(j['seq']))
			ds = (j['dotstr'].count('(')+j['dotstr'].count(')')) / float(len(j['seq']))
			v_auc_ls = []
			gini_ls = []
			for n in xrange(100):
				v_ls = []
				for d in list(dot):
					if d in ['.']:
						v = random.uniform(ss_min, ss_max)
					elif d in ['(', ')']:
						v = random.uniform(ds_min, ds_max)
					else:
						v = random.uniform(0, 1)
					v_ls.append(v)
				ct_list = map(int, list(dot.replace('.','1').replace('(','0').replace(')','0')))
				fpr, tpr, _ = roc_curve(ct_list, v_ls)
				gini = gj.gini(v_ls)
				v_auc = auc(fpr, tpr)
				v_auc_ls.append(v_auc)
				gini_ls.append(gini)
			print>>SAVEFN, '\t'.join(map(str, [i, ss,ds,','.join(map(str, v_auc_ls)), ','.join(map(str, gini_ls)),dot, ','.join(map(str, v_ls))]))



dot = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/mouse.dot'
savefn = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/simulate_shape/mouse.simulate.shape.txt'
simulate_dot(dot=dot, savefn=savefn)
