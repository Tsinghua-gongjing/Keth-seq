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

def read_icshape_out(out=None, pureID=0):
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
			out_dict[tx_id]['length'] = length
			out_dict[tx_id]['rpkm'] = rpkm
			out_dict[tx_id]['reactivity_ls'] = reactivity_ls
	gj.printFuncRun('read_icshape_out')
	return out_dict.to_dict()

def calc_auc(file_list, df_list, ct_list, access_ls, seq, savefn, accessibility_cutoff=3, reactivity_cutoff=0, base_used='G', plot_ss_ds=0):
    file_list = file_list[0:]
    df_list = df_list[0:]
    
    fpr_list = []
    tpr_list = []
    roc_auc_list = []
    
    for i in xrange(len(df_list)):
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
    plt.title('base: %s, access: %s, reactivity: %s'%(base_used, accessibility_cutoff, reactivity_cutoff))
    plt.legend(loc=4)
    if not savefn is None:
        plt.savefig(savefn)
    else:
        plt.show()

    return fpr_list, tpr_list, roc_auc_list, select

def read_score(score_tab, ref):
    fa = Fasta(ref)
    fa_dict = {}
    for i,j in fa.items():
        fa_dict[i.split('\t')[0]] = j
    
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

def shape_structure(shape, dot, savedir, base='G'):
	if shape.endswith('.out'):
		out_dict = read_icshape_out(shape)
	elif shape.endswith('score.tab'):
		out_dict = read_score(shape, '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/human.dot.dedup.fa')
	seq_dict = read_dot(dot)
	print out_dict.keys()

	seq_id_ls = []
	roc_auc_ls = []
	base_num_ls = []

	for i,j in out_dict.items():
		if seq_dict.has_key(i):
			print "calc_auc: %s"%(i)

			reactivity_ls = [np.nan if v == 'NULL' else float(v) for v in j['reactivity_ls']]
			reactivity_df = pd.DataFrame({'reactivity':reactivity_ls})
			ct_list = list(seq_dict[i]['dotstr'].replace('.','1').replace('(','0').replace(')','0'))
			print reactivity_ls
			if reactivity_ls.count(np.nan) == len(reactivity_ls):
				continue
			reactivity_ls_no_null = [r for r in reactivity_ls if not np.isnan(r)]
			if len(reactivity_ls_no_null) <= 5:
				continue
			reactivity_ls_base_no_null = [r for r,s in zip(reactivity_ls, list(seq_dict[i]['seq'])) if (not np.isnan(r)) and (s in list(base))]
			if len(reactivity_ls_base_no_null) <= 2:
				continue

			print "df", reactivity_df.shape
			print "dotstr", len(ct_list)
			print "seq", len(seq_dict[i]['seq']), list(seq_dict[i]['seq'])
			savefn = '%s/%s.png'%(savedir, i.split('/')[0])
			fpr_list, tpr_list, roc_auc_list, select = calc_auc(file_list=[i], 
					 df_list=[reactivity_df], 
					 ct_list=ct_list, 
					 access_ls=None, 
					 seq=seq_dict[i]['seq'], 
					 savefn=savefn, 
					 accessibility_cutoff=None, 
					 reactivity_cutoff=0, 
					 base_used=base, 
					 plot_ss_ds=0)
			print "auc", roc_auc_list
			roc_auc_ls.append(roc_auc_list[0])
			seq_id_ls.append(i)
			base_num_ls.append(select.shape[0])
	roc_auc_df = pd.DataFrame({'seq_id':seq_id_ls, 'AUC':roc_auc_ls, 'base_num':base_num_ls})
	roc_auc_df.to_csv(savedir+'/AUC.txt', header=True, index=False, sep='\t')
	return roc_auc_ls

shape = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.rfam.T1t10.out'
dot = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/mouse.dot'
savedir = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/mouse'
# shape_structure(shape, dot, savedir, base='G')

shape = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mes/wc/out/shape.out'
dot = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/mouse.dot'
savedir = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mes/wc/out/mouse'
# shape_structure(shape, dot, savedir, base='ATCG')

shape = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_kethoxal.rfam.T1t200.out'
dot = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/human.dot'
savedir = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/keth_hela'
# shape_structure(shape, dot, savedir, base='G')

shape = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/DMS_Nature_2014_Weissman/raw/merge/fibroblast_vivo.rfam.score.tab'
dot = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/human.dot'
savedir = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/DMSseq_fibroblast_vivo'
# shape_structure(shape, dot, savedir, base='AC')

shape = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/DMS_Nature_2014_Weissman/raw/merge/K562_vivo.rfam.score.tab'
dot = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/human.dot'
savedir = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/DMSseq_K562_vivo'
# shape_structure(shape, dot, savedir, base='AC')

def compare_AUC(AUC1, AUC2, method1='keth-seq', method2='icSHAPE', savefn='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/keth_vs_icshape.pdf'):
	df_auc1 = pd.read_csv(AUC1, header=0, sep='\t')
	df_auc1.columns = ['AUC(%s)'%(method1), 'base_num(%s)'%(method1), 'seq_id']
	df_auc2 = pd.read_csv(AUC2, header=0, sep='\t')
	df_auc2.columns = ['AUC(%s)'%(method2), 'base_num(%s)'%(method2), 'seq_id']
	print df_auc1
	print df_auc2
	df_AUC = df_auc1.merge(df_auc2, on='seq_id')
	print df_AUC
	df_AUC_common = df_AUC.dropna(how='any')
	print df_AUC_common

	fig,ax=plt.subplots(figsize=(8,8))
	sns.scatterplot(x='AUC(%s)'%(method1), y='AUC(%s)'%(method2), data=df_AUC_common)
	plt.plot([-0.05, 1.05], [-0.05, 1.05], color='black', linestyle='--')
	ax.axis('square')
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.tight_layout()
	# savefn = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/keth_vs_icshape.pdf'
	plt.savefig(savefn)

	# df_AUC_common.to_csv(savefn.replace('.pdf', '.txt'), header=True, index=False, sep='\t')


# AUC1 = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/mouse/AUC.txt'
# AUC2 = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mes/wc/out/mouse/AUC.txt'
# compare_AUC(AUC1, AUC2)


AUC1 = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/F2c.mes_kethseq.AUC.txt'
AUC2 = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/F2c.mes_icSHAPE.AUC.txt'
compare_AUC(AUC1, AUC2, savefn='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/Keth-seq/results/F2c.keth_vs_icshape.AUC.pdf')



