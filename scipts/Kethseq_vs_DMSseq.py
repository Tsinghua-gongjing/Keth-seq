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
from sklearn.metrics import roc_curve, auc

def calc_auc2(file_list, df_list, ct_list, access_ls, seq, savefn, accessibility_cutoff=3, reactivity_cutoff=0, base_used_ls=None, plot_ss_ds=0):
    file_list = file_list[0:]
    df_list = df_list[0:]
    
    fpr_list = []
    tpr_list = []
    roc_auc_list = []
    
    for i,base_used in zip(xrange(len(df_list)), base_used_ls):
        df_list[i]['ss_ds'] = map(int, ct_list)
        df_list[i]['base'] = list(seq)
        if not access_ls is None:
            df_list[i]['accessibility'] = access_ls
        clean = df_list[i].dropna()
        if not accessibility_cutoff is None:
            select = clean[clean.accessibility >= accessibility_cutoff]
        else:
            select = clean.loc[:,:]
        select = select[select.reactivity >= reactivity_cutoff]
        select = select[select.base.isin(list(base_used))]
        
        if plot_ss_ds:
            single_strand = select[select.ss_ds == 1]
            double_strand = select[select.ss_ds == 0]
            print select.head(5)
            print select['ss_ds'].value_counts()
            print single_strand.describe()
            print double_strand.describe()
            sns.distplot(single_strand.dropna()['reactivity'], norm_hist=True, kde=False, bins=10)
            sns.distplot(double_strand.dropna()['reactivity'], norm_hist=True, kde=False, bins=10)
            plt.title(file_list[i])
            plt.show()
            
        fpr, tpr, _ = roc_curve(select['ss_ds'], select['reactivity'])
        fpr_list.append(fpr)
        tpr_list.append(tpr)
        roc_auc_list.append(auc(fpr, tpr))
        
    plt.close()
    plt.figure()
    lw = 2
    with sns.axes_style("ticks"):
        fig,ax=plt.subplots(figsize=(8,8))
        for i in xrange(len(df_list)):
            print "plot: %s, %s"%(i, file_list[i])
            if i > 10:
                continue
            plt.plot(fpr_list[i],tpr_list[i],lw = lw,label = '%s (area = %0.2f)' %(file_list[i],roc_auc_list[i]))
            
    plt.plot([-0.05, 1.05], [-0.05, 1.05], color='black', lw=lw, linestyle='--',label='Luck')
    ax.axis('square')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('access: %s, reactivity: %s'%(accessibility_cutoff, reactivity_cutoff))
    plt.legend(loc=4)
#     plt.legend(bbox_to_anchor=(1, 1), loc=2)
    if not savefn is None:
        plt.savefig(savefn)
    else:
        plt.show()

    return fpr_list, tpr_list, roc_auc_list


def read_score(score_tab, ref):
    fa = Fasta(ref)
    fa_dict = {}
    for i,j in fa.items():
        fa_dict[i.split(' ')[0]] = j
    
    score_dict = nested_dict()
    with open(score_tab, 'r') as TXT:
        for line in TXT:
            line = line.strip()
            if not line or line.startswith('@'): continue
            arr = line.split('\t')
            if arr[1] == '-': continue
            score_dict[arr[0]][int(arr[2])] = arr[7]
    score_dict = score_dict.to_dict()
#     print score_dict
    
    reactivity_dict = nested_dict(2, list)
    for i,j in score_dict.items():
        for p in xrange(1, len(fa_dict[i])+1):
            if p not in score_dict[i]:
                r = 'NULL'
            elif score_dict[i][p] == '-1':
                r = 'NULL'
            else:
                r = score_dict[i][p]
            reactivity_dict[i]['reactivity_ls'].append(r)
            
    return reactivity_dict.to_dict()

fibroblast_vivo_reactivity_dict = read_score('./DMSseq/fibroblast_vivo.rRNA.score.tab', './human_ribosomalRNA_4.fa')
fibroblast_vitro_reactivity_dict = read_score('./DMSseq/fibroblast_vitro.rRNA.score.tab', './human_ribosomalRNA_4.fa')
K562_vivo_reactivity_dict = read_score('./DMSseq/K562_vivo.rRNA.score.tab', './human_ribosomalRNA_4.fa')
K562_vitro_reactivity_dict = read_score('./DMSseq/K562_vitro.rRNA.score.tab', './human_ribosomalRNA_4.fa')


fibroblast_vivo_18S = fibroblast_vivo_reactivity_dict['18S']['reactivity_ls']
fibroblast_vivo_18S = [np.nan if i == 'NULL' else float(i) for i in fibroblast_vivo_18S]
fibroblast_vivo_18S_df = pd.DataFrame({'reactivity':fibroblast_vivo_18S})

fibroblast_vitro_18S = fibroblast_vitro_reactivity_dict['18S']['reactivity_ls']
fibroblast_vitro_18S = [np.nan if i == 'NULL' else float(i) for i in fibroblast_vitro_18S]
fibroblast_vitro_18S_df = pd.DataFrame({'reactivity':fibroblast_vitro_18S})

K562_vivo_18S = K562_vivo_reactivity_dict['18S']['reactivity_ls']
K562_vivo_18S = [np.nan if i == 'NULL' else float(i) for i in K562_vivo_18S]
K562_vivo_18S_df = pd.DataFrame({'reactivity':K562_vivo_18S})

K562_vitro_18S = K562_vitro_reactivity_dict['18S']['reactivity_ls']
K562_vitro_18S = [np.nan if i == 'NULL' else float(i) for i in K562_vitro_18S]
K562_vitro_18S_df = pd.DataFrame({'reactivity':K562_vitro_18S})


file_list = ['keth-seq(HeLa)', 'DMSseq(F,vivo)',  'DMSseq(K562,vivo)']
df_list = [Hela_kethoxal_18S_df,  
           fibroblast_vivo_18S_df, K562_vivo_18S_df, ]
base_used_ls = ['G','AC', 'AC']

accessibility_cutoff_ls = [3]#[0, 1, 3, 5, 10, 20, 30]
reactivity_cutoff_ls = [0]#[0, 0.01, 0.05]

for accessibility_cutoff in accessibility_cutoff_ls:
    for reactivity_cutoff in reactivity_cutoff_ls:
            savefn_dir = '/Users/gongjing/SeafileSyn/github/KS_test/scripts/AUC/compare_auc_plots'
            savefn = '%s/human.18S.base%s.access%s.reactivity%s.pdf'%(savefn_dir, 'Sep',accessibility_cutoff, reactivity_cutoff, )
            fpr_list, tpr_list, roc_auc_list = calc_auc2(file_list, df_list, dot_dict_human18S_new, list(accessibility['accessibility']), 
                    fa_dict_rRNA_human['18S'], savefn, accessibility_cutoff, reactivity_cutoff, base_used_ls, 0)