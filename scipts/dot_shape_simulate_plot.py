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

df_simu = pd.read_csv('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/SF10c.mouse.simulate.shape.txt', header=None, sep='\t')
df_simu.columns = ['id', 'ratio(ss)', 'ratio(ds)', 'auc_ls', 'gini_ls', 'dot', 'reactivity']
df_simu.sort_values(by='ratio(ss)', inplace=True)

df_simu.head()

df_simu['gini(mean)'] = [np.mean(map(float,i.split(','))) for i in df_simu['gini_ls']]

# fig,ax=plt.subplots(figsize=(16,5))
# df_t = df_simu.loc[:,['ratio(ss)', 'ratio(ds)']]
# print df_t.shape
# df_t.plot.bar(stacked=True, ax=ax, linewidth=0, width=1)
# ax.set_xticks([])
# plt.legend(bbox_to_anchor=(1, 1), loc=2)
# ax.set_ylim(0,1)
# ax.set_ylabel('Percentage')
# ax.set_xlabel('RNA with known structure (from Rfam)')

# plt.tight_layout()
# plt.savefig('./mouse.ss_ds.pdf')

fig,ax=plt.subplots()
sns.jointplot(x='ratio(ds)',y='gini(mean)', data=df_simu, stat_func=stats.pearsonr,kind="reg", size=7)
plt.tight_layout()
plt.savefig('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/SF10c.mouse.simulate.shape.dsRatio_vs_gini.pdf')