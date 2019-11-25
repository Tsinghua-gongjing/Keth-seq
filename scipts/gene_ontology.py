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


BP_vitro = pd.read_csv('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/SF12c.in_vitro_95_GO_BP.txt', header=0, sep='\t')
BP_vitro = BP_vitro[BP_vitro['PValue']<0.05]
BP_vitro['-log10(PValue)'] = -np.log10(BP_vitro['PValue'])
BP_vitro

fig,ax=plt.subplots(figsize=(20, 6))

sns.barplot(y='Term', x='-log10(PValue)', data=BP_vitro.iloc[0:10,:], color='#8172B2')
plt.tight_layout()
plt.savefig('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/SF12c.in_vitro_95_GO_BP.pdf')
plt.close()


BP_vivo = pd.read_csv('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/SF12f.in_vivo_105_GO_BP.txt', header=0, sep='\t')
BP_vivo = BP_vivo[BP_vivo['PValue']<0.005]
BP_vivo['-log10(PValue)'] = -np.log10(BP_vivo['PValue'])
BP_vivo

fig,ax=plt.subplots(figsize=(20, 6))

sns.barplot(y='Term', x='-log10(PValue)', data=BP_vivo.iloc[0:10,:], color='#8172B2')
plt.tight_layout()
plt.savefig('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/SF12f.in_vivo_105_GO_BP.pdf')
plt.close()