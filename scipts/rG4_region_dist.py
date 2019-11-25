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

fig, ax = plt.subplots(1,2,figsize=(8,4))

# in vitro, our control
sizes = [11, 34, 35, 15]
labels = ["UTR'5", 'CDS', "UTR'3", 'Others']
total = sum(sizes)
ax[0].pie(sizes, labels=labels, 
        autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)
ax[0].axis('equal')


# in vivo
sizes = [12, 52, 34, 7]
labels = ["UTR'5", 'CDS', "UTR'3", 'Others']
total = sum(sizes)
ax[1].pie(sizes, labels=labels, 
        autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)
ax[1].axis('equal')

# in vitro, Rg4-seq control
# sizes = [14, 56, 62, 15]
# labels = ["UTR'5", 'CDS', "UTR'3", 'Others']
# total = sum(sizes)
# ax[2].pie(sizes, labels=labels, 
#         autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)
# ax[2].axis('equal')

plt.savefig('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/Sf12be.in_vivo_vitro_rG4_element_dist_pie.pdf')
plt.show()