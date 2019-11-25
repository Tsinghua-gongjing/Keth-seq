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


def read_icshape_out(out=None, pureID=1):
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
            out_dict[tx_id]['length'] = length
            out_dict[tx_id]['rpkm'] = rpkm
            out_dict[tx_id]['reactivity_ls'] = reactivity_ls
    return out_dict

def read_dot(dot=None):
    dot_dict = nested_dict()
    with open(dot, 'r') as DOT:
        for n,line in enumerate(DOT):
            line = line.strip()
            if not line or line.startswith('#'): continue
            if n == 0:
                dot_dict['name'] = line.replace('>','')
            if n == 1:
                dot_dict['seq'] = line
            if n == 2:
                dot_dict['dotbracket'] = line
    #print dot_dict
    return dot_dict

def read_mouse_dot():
    dot_dict_18S = read_dot('../data/rRNA/mouse_RNA_strand_16S.dot')
    mouse_dot_18S = list(dot_dict_18S['dotbracket'])
    mouse_dot_18S_new = ''.join(mouse_dot_18S[1:][0:18] + ['.'] + mouse_dot_18S[1:][18:] + ['.'])
    print mouse_dot_18S_new, len(mouse_dot_18S_new)

    mouse_ct_18S = list(dot_dict_18S['dotbracket'].replace('.','1').replace('(','0').replace(')','0'))

    print "original", len(mouse_ct_18S)
    # convert 1869bp to 1870bp ct dot file
    # this is based on sequence alignment
    mouse_ct_18S_new = map(int, mouse_ct_18S[1:][0:18] + ['1'] + mouse_ct_18S[1:][18:] + ['1'])
    print "converted", len(mouse_ct_18S_new)


    return mouse_dot_18S_new, mouse_ct_18S_new

def read_access():
    accessibility = pd.read_csv('../data/rRNA/18s_o2_sasa.sort.txt', header=None, sep='\t')
    accessibility.columns = ['position', 'base', 'accessibility']
    accessibility.head()

    # seq alignment: mouse(1870), human(1869)
    m = 'ACCTGGTTGATCCTGCCAGGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCCTTCCCGGGGGGGGATGCGTGCATTTATCAGATCAAAACCAACCCGGTGAGCTCCCTC-CCGGCTCCGGCCGGGGGTCGGGCGCCGGCGGC-TTGGTGACTCTAGATAACCTCGGGCCGATCGCACGCCCCCCGTGGCGGCGACGACCCATTCGAACGTCTGCCCTATCAACTTTCGATGGTAGTCGCCGTGCCTACCATGGTGACCACGGGTGACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCACTCCCGACCCGGGGAGGTAGTGACGAAAAATAACAATACAGGACTCTTTCGAGGCCCTGTAATTGGAATGAGTCCACTTTAAATCCTTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGTCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACTGGCGATGCGGCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGG-GCCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGATTCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTGCGCCGGCAGGCGCGGGTAACCCGTTGAACCCCATTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGAATTCCCAGTAAGTGCGGGTCATAAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGCCCTCGGATCGGCCCCGCCGGGGTCGGCCCACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAACTTGACTATCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA'
    h = 'ACCTGGTTGATCCTGCCA-GTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTGAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCC-TTCGC-GGGGGGGATGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGC-CCCTCTCCGGCCCCGGCCGGGGGGCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGCCCCCCGTGGCGGCGACGACCCATTCGAACGTCTGCCCTATCAACTTTCGATGGTAGTCGCCGTGCCTACCATGGTGACCACGGGTGACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCACTCCCGACCCGGGGAGGTAGTGACGAAAAATAACAATACAGGACTCTTTCGAGGCCCTGTAATTGGAATGAGTCCACTTTAAATCCTTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGCCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACCGGCGATGCGGCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGATTCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACGCCGGCAGGCGCGGGTAACCCGTTGAACCCCATTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGAATTCCCAGTAAGTGCGGGTCATAAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGCCCTCGGATCGGCCCCGCCGGGGTCGGCCCACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAACTTGACTATCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA'

    a_ls = list(accessibility['accessibility'])
    print "original", len(a_ls)
    new_ls = []

    n = 0
    for i,j in zip(m, h):
        if j != '-':
            n += 1
        if i != '-':
            if j == '-':
                new_ls.append('NULL')
            else:
                new_ls.append(a_ls[n])
        
    print "converted 1", len(new_ls)

    # convert human 18S 1869bp to mouse 18S 1870bp (accessibility info)
    mouse_acces_ls = new_ls + ['NULL']
    mouse_acces_ls = [np.nan if i == 'NULL' else i for i in mouse_acces_ls]
    print "converted 2", len(mouse_acces_ls)

    return mouse_acces_ls

out_dict_kethoxal_rRNA_noTreat = read_icshape_out('../data/shape/mES_kethoxal_vs_no-treat.rRNA.shape.out.txt')
kethoxal_notreat = out_dict_kethoxal_rRNA_noTreat['18S']['reactivity_ls']
kethoxal_notreat = [np.nan if i == 'NULL' else float(i) for i in kethoxal_notreat]
kethoxal_notreat_df = pd.DataFrame({'reactivity':kethoxal_notreat})


out_dict_icshape_mES_panpan_rRNA = read_icshape_out('../data/icSHAPE/mes.rRNA.out')
# based on seq alignment, convert signal from 1869 to 1870
icshape_panpan_1869 = out_dict_icshape_mES_panpan_rRNA['mouse_small']['reactivity_ls']
icshape_panpan_1870 = icshape_panpan_1869[1:][0:18] + ['NULL'] + icshape_panpan_1869[1:][18:] + ['NULL']
icshape_panpan_1870 = [np.nan if i == 'NULL' else float(i) for i in icshape_panpan_1870]
fromPan1870_df = pd.DataFrame({'reactivity':icshape_panpan_1870})


mouse_acces_ls = read_access()
mouse_dot_18S_new, mouse_ct_18S_new = read_mouse_dot()


fa_dict1 = Fasta('../data/rRNA/mouse_ribosomalRNA_4.fa', key_fn=lambda key:key.split("\t")[0])
fa_dict_rRNA_mouse = {i.split()[0].split('.')[0]:j[0:] for i,j in fa_dict1.items()}
# fa_dict_rRNA_mouse.keys()


df_reactivity_df= pd.DataFrame({'base':list(fa_dict_rRNA_mouse['18S']),
                               'dot':list(mouse_dot_18S_new),
                               'keth-seq':kethoxal_notreat_df['reactivity'],
                               'icSHAPE':fromPan1870_df['reactivity']})
df_reactivity_df.head()


### plot diff of single stranded bases
df_reactivity_df = df_reactivity_df[df_reactivity_df['base']=='G']
df_reactivity_df = df_reactivity_df[df_reactivity_df['dot'].isin(['.', ])]
df_reactivity_df.dropna(how='any', inplace=True)
df_reactivity_df['diff'] = df_reactivity_df['keth-seq'] - df_reactivity_df['icSHAPE']
df_reactivity_df['pos'] = range(df_reactivity_df.shape[0])
df_reactivity_df

fig,ax=plt.subplots(2,1, figsize=(20, 12), sharex=True)
df_reactivity_df.plot(x='pos', y='keth-seq', ax=ax[0], label='keth-seq')
df_reactivity_df.plot(x='pos', y='icSHAPE', ax=ax[0], label='icSHAPE')
df_reactivity_df.plot(x='pos', y='diff', ax=ax[1], label='diff', color='green')
ax[1].set_ylim(-1.05, 1.05)
plt.axhline(y=0, xmin=0, xmax=1, hold=None, color='lightgrey', ls='--')
ax[0].set_ylabel('Reactivity')
ax[1].set_ylabel('Reactivity difference \n(keth-seq - icSHAPE)')
plt.tight_layout()
savefn = '../results/SF10a.18S_G_base_reactivity_diff.pdf'
plt.savefig(savefn)
plt.close()


### plot AUC
file_list = ['kethoxal', 'icshape2_panpan']
df_list = [kethoxal_notreat_df, fromPan1870_df]

accessibility_cutoff_ls = [3,]
reactivity_cutoff_ls = [0, ]
base_used_ls = ['G', 'ATCG']

for accessibility_cutoff in accessibility_cutoff_ls:
    for reactivity_cutoff in reactivity_cutoff_ls:
        for base_used in base_used_ls:
            savefn_dir = '../results'
            savefn = '%s/SF10a.mouse.18S.base%s.access%s.reactivity%s.pdf'%(savefn_dir, 'Sep',accessibility_cutoff, reactivity_cutoff, )
            fpr_list, tpr_list, roc_auc_list = calc_auc2(file_list, df_list, mouse_ct_18S_new, mouse_acces_ls, 
                    fa_dict_rRNA_mouse['18S'], savefn, accessibility_cutoff, reactivity_cutoff, base_used_ls, 0)
