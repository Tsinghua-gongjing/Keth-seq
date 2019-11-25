import gj
import argparse
import subprocess
from nested_dict import nested_dict
import pandas as pd 
import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_context("poster")
sns.set_style("ticks")

def read_collapse(fastq=None):
	gj.printFuncRun('read_collapse')
	gj.printFuncArgs()
	collapse_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/readCollapse.pl'
	collapse_fq = fastq.replace('fastq','rmdup.fastq')
	seq_freq_fa = fastq.replace('fastq','fa')
	subprocess.call(["%s -U %s -o %s -f %s"%(collapse_pl, fastq, collapse_fq, seq_freq_fa)],shell=True)
	gj.printFuncRun('read_collapse')
	return collapse_fq

def read_collapse_PE(fastq1=None, fastq2=None):
	gj.printFuncRun('read_collapse_PE')
	gj.printFuncArgs()
	collapse_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/readCollapse.pl'
	collapse_fq1 = fastq1.replace('fastq','rmdup.fastq')
	seq_freq_fa1 = fastq1.replace('fastq','fa')
	collapse_fq2 = fastq2.replace('fastq','rmdup.fastq')
	seq_freq_fa2 = fastq2.replace('fastq','fa')
	subprocess.call(["%s -U %s -o %s -f %s"%(collapse_pl, fastq1, collapse_fq1, seq_freq_fa1)],shell=True)
	subprocess.call(["%s -U %s -o %s -f %s"%(collapse_pl, fastq2, collapse_fq2, seq_freq_fa2)],shell=True)
	gj.printFuncRun('read_collapse_PE')

def remove_adapter(fastq=None, trimmed_fastq=None):
	gj.printFuncRun('remove_adapter')
	gj.printFuncArgs()
	trimming_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/trimming.pl'
	if trimmed_fastq is None:
		trimmed_fastq = fastq.replace('rmdup.fastq','trimmed.fastq')
	# adapter_fa = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/data/adapter/kethoxal.fa'
	# subprocess.call(["%s -U %s -o %s -l 13 -t 0 -c phred33 -a %s -m 0"%(trimming_pl, fastq, trimmed_fastq, adapter_fa)],shell=True)

	# for lib construct as Yichengqi's lab
	adapter_fa = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/data/adapter/kethoxal_yichengqi.fa'
	subprocess.call(["%s -U %s -o %s -l 10 -t 0 -c phred33 -a %s -m 0"%(trimming_pl, fastq, trimmed_fastq, adapter_fa)],shell=True)

	"""
	for min_len in [50,25]:
		trimmed_fastq_minLen = trimmed_fastq.replace('fastq','minLen'+str(min_len)+'.fastq')
		subprocess.call(["%s -U %s -o %s -l 13 -t 0 -c phred33 -a %s -m %s"%(trimming_pl, fastq, trimmed_fastq, adapter_fa, min_len)],shell=True)
	"""
	gj.printFuncRun('remove_adapter')
	return trimmed_fastq

def remove_adapter_PE(fastq1=None, fastq2=None):
	gj.printFuncRun('remove_adapter_PE')
	gj.printFuncArgs()
	trimming_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/trimming.pl'
	trimmed_fastq1 = fastq1.replace('rmdup.fastq','trimmed.fastq')
	trimmed_fastq2 = fastq2.replace('rmdup.fastq','trimmed.fastq')
	adapter_fa = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/N02381_CY_80-65347437_DNASEQ/PF_data/170705-X1B/TruSeq3-PE.fa'
	subprocess.call(["%s -U %s -o %s -l 13 -t 0 -c phred33 -a %s -m 50"%(trimming_pl, fastq1, trimmed_fastq1, adapter_fa)],shell=True)
	subprocess.call(["%s -U %s -o %s -l 13 -t 0 -c phred33 -a %s -m 50"%(trimming_pl, fastq2, trimmed_fastq2, adapter_fa)],shell=True)
	gj.printFuncRun('remove_adapter_PE')

def remove_adapter_PE_new(fastq1=None, fastq2=None):
	gj.printFuncRun('remove_adapter_PE_new')
	gj.printFuncArgs()
	trimmed_fastq1 = fastq1.replace('.fastq','trimmed.fastq')
	trimmed_fastq2 = fastq2.replace('.fastq','trimmed.fastq')
	adapter_fa1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/data/adapter/kethoxal.fa'
	adapter_fa2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/data/adapter/kethoxal_rev.fa'
	adapter_fa = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/data/adapter/kethoxal_PE.fa'
	trimmomatic = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/bin/trimmomatic-0.30.jar'
	R1_paired = fastq1.replace('fastq', 'paired.fastq')
	R1_unpaired = fastq1.replace('fastq', 'unpaired.fastq')
	R2_paired = fastq2.replace('fastq', 'paired.fastq')
	R2_unpaired = fastq2.replace('fastq', 'unpaired.fastq')
	trimlog = fastq1.replace('_R1_001.fastq', '.trim.log')
	subprocess.call(["java -jar %s PE -threads 32 -phred33 -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:33"%(trimmomatic, trimlog, fastq1, fastq2, R1_paired, R1_unpaired, R2_paired, R2_unpaired, adapter_fa)],shell=True)

	clip_fastq1 = R1_paired.replace('fastq', 'clip.fastq')
	clip_trimlog1 = clip_fastq1 + '.log'
	clip_fastq2 = R2_paired.replace('fastq', 'clip.fastq')
	clip_trimlog2 = clip_fastq2 + '.log'
	subprocess.call(["cutadapt -u 13 -o %s %s"%(clip_fastq1, trimmed_fastq1)], shell=True)
	subprocess.call(["cutadapt -u -13 -o %s %s"%(clip_fastq2, trimmed_fastq2)], shell=True)

	gj.printFuncRun('remove_adapter_PE_new')

def mapping(fastq=None, mapper='bowtie2', index_dir=None):
	gj.printFuncRun('mapping')
	gj.printFuncArgs()
	map_sam = fastq.replace('trimmed.fastq','sam')
	map_rRNA_sam = fastq.replace('trimmed.fastq','rRNA.sam')

	#index = '/Share/home/zhangqf/database/GenomeAnnotation/INDEX/Bowtie2/mm_rRNA/mm_rRNA'
	# index = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/rRNA/index/rrna_uniq'
	index = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/hg38/transcriptome/ribosomalRNAindex/ribosomalRNA'
	rRNA_unmap_fastq = fastq.replace('fastq','rRNAUnmap.fastq')
	subprocess.call(["bowtie2 -U %s -S %s -x %s --non-deterministic --time --un %s"%(fastq, map_rRNA_sam, index, rRNA_unmap_fastq)],shell=True)
	
	index = '/Share/home/zhangqf/database/GenomeAnnotation/INDEX/Bowtie2/hg38/Gencode_transcriptome/whole_transcriptome/hg38'
	subprocess.call(["bowtie2 -U %s -S %s -x %s --non-deterministic --time"%(rRNA_unmap_fastq, map_sam, index)],shell=True)

	gj.printFuncRun('mapping')
	return map_sam

def map_rfam(fastq=None, map_sam=None, index=None):
	gj.printFuncRun('mapping')
	gj.printFuncArgs()
	subprocess.call(["bowtie2 -U %s -S %s -x %s --non-deterministic --time"%(fastq, map_sam, index)],shell=True)
	gj.printFuncRun('mapping')
	return map_sam

def mapping_PE(fastq1=None, fastq2=None, mapper='bowtie2', index_dir=None):
	gj.printFuncRun('mapping_PE')
	gj.printFuncArgs()
	map_sam = fastq1.replace('paired.clip.fastq','sam')
	map_rRNA_sam = fastq1.replace('paired.clip.fastq','rRNA.sam')

	# map to rRNA
	index = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/rRNA/index/rrna_uniq'
	rRNA_unmap_fastq = fastq1.replace('fastq','rRNAUnmap.fastq')
	subprocess.call(["bowtie2 -p 12 -1 %s -2 %s -x %s -S %s --non-deterministic --time --un-conc %s --no-unal"%(fastq1, fastq2, index, map_rRNA_sam, rRNA_unmap_fastq)],shell=True)


	# map to transcriptome
	index = '/Share/home/zhangqf/database/GenomeAnnotation/INDEX/Bowtie2/mm10/Gencode_transcriptome/whole_transcriptome/mm10'
	#index = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Cyprinus_carpio/rna'
	#map_sam = fastq1.replace()
	rRNA_unmap_fastq1 = fastq1.replace('fastq','rRNAUnmap.1.fastq')
	rRNA_unmap_fastq2 = fastq1.replace('fastq','rRNAUnmap.2.fastq')
	subprocess.call(["bowtie2 -p 12 -1 %s -2 %s -x %s -S %s --non-deterministic --time --no-unal"%(rRNA_unmap_fastq1, rRNA_unmap_fastq2, index, map_sam)],shell=True)

	gj.printFuncRun('mapping_PE')

def rpkm_cal(sam=None):
	gj.printFuncRun('rpkm_cal')
	gj.printFuncArgs()
	estimateRPKM_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/estimateRPKM.pl'
	sam_rpkm = sam.replace('sam','rpkm')
	subprocess.call(["%s -i %s -o %s"%(estimateRPKM_pl, sam, sam_rpkm)],shell=True)
	gj.printFuncRun('rpkm_cal')
	return sam_rpkm

def RT_cal(sam=None):
	gj.printFuncRun('RT_cal')
	gj.printFuncArgs()
	sam_rt = sam.replace('sam', 'rt')
	sam_rpkm = sam.replace('sam','rpkm')
	calcRT_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/calcRT.pl'
	subprocess.call(["%s -i %s -o %s -r %s -c 1"%(calcRT_pl, sam, sam_rt, sam_rpkm)],shell=True)
	gj.printFuncRun('RT_cal')
	return sam_rt

def read_clean_map_rt(fastq):
	gj.printFuncRun('read_clean_map_rt')
	gj.printFuncArgs()
	collapse_fq = read_collapse(fastq=fastq)
	trimmed_fastq = remove_adapter(fastq=collapse_fq)
	map_sam = mapping(fastq=trimmed_fastq)
	sam_rpkm = rpkm_cal(sam=map_sam)
	sam_rt = RT_cal(sam=map_sam)
	gj.printFuncRun('read_clean_map_rt')

def RT_correlation(rt1=None, rt2=None, rt_corr=None, coverage_cutoff=0, background_base_density=0):
	gj.printFuncRun('RT_correlation')
	gj.printFuncArgs()
	correlationRT_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/correlationRT.pl'
	subprocess.call(["perl %s -1 %s -2 %s -T %s -b %s > %s"%(correlationRT_pl, rt1, rt2, coverage_cutoff, background_base_density, rt_corr)],shell=True)
	gj.printFuncRun('RT_correlation')

def RT_combine(rt1=None, rt2=None, rt_comb=None):
	gj.printFuncRun('RT_combine')
	gj.printFuncArgs()
	combineRT_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/combineRTreplicates.pl'
	subprocess.call(["perl %s -i %s:%s -o %s"%(combineRT_pl, rt1, rt2, rt_comb)],shell=True)
	gj.printFuncRun('RT_combine')

def RT_normalize(rt=None):
	gj.printFuncRun('RT_normalize')
	gj.printFuncArgs()
	normalize_rt = rt.replace('.rt', '.normalized.rt').replace('RT', 'normalized.RT')
	normalizeRT_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/normalizeRTfile.pl'
	subprocess.call(["perl %s -i %s -o %s -m mean:vigintile2 -d 32 -l 32"%(normalizeRT_pl, rt, normalize_rt)],shell=True)
	gj.printFuncRun('RT_normalize')

def calc_enrich(f_normalized_rt=None, b_normalized_rt=None, icshape_tmp_out=None, x=0.25):
	gj.printFuncRun('calc_enrich')
	gj.printFuncArgs()
	calc_enrich_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/calcEnrich.pl'
	subprocess.call(["perl %s -f %s -b %s -o %s -w factor5:scaling1 -x %s"%(calc_enrich_pl, f_normalized_rt, b_normalized_rt, icshape_tmp_out, x)],shell=True)
	gj.printFuncRun('calc_enrich')

def filter_enrich(icshape_tmp_out=None, average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30):
	gj.printFuncRun('filter_enrich')
	gj.printFuncArgs()
	icshape_out = icshape_tmp_out.replace('.tmp', '.T%st%s'%(average_coverage, background_base_density))
	filter_enrich_pl = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/scripts/icSHAPE-master/scripts/filterEnrich.pl'
	subprocess.call(["perl %s -i %s -o %s -T %s -t %s -s %s -e %s"%(filter_enrich_pl, icshape_tmp_out, icshape_out, average_coverage, background_base_density, skip_leading, skilp_tailing)],shell=True)
	gj.printFuncRun('filter_enrich')

def library_info2():
	library_info_dict = nested_dict()
	library_info_dict['dir_prefix'] = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data'
	library_info_dict['lib1']['dir'] = '16-03-25_library_of_concentration'
	library_info_dict['lib2']['dir'] = '16-08-08_16_library_invivo_invitro'
	library_info_dict['lib3']['dir'] = '16-09-09_10_library_G4-glucose_time'
	library_info_dict['lib4']['dir'] = '16-10-14_4_library_08-08_rerun'
	library_info_dict['lib5']['dir'] = '16-11-14_7_library_total_Kethoxal_remove'
	library_info_dict['lib6']['dir'] = '19-04-03'

	library_info_dict['lib5']['fq'] = {'control1':'CHe-XC-M1C_S1_L006_R1_001.fastq', 'control2':'CHe-XC-M2C_S2_L006_R1_001.fastq', 
								 'kethoxal1':'CHe-XC-M1K_S3_L006_R1_001.fastq', 'kethoxal2':'CHe-XC-M2K_S4_L006_R1_001.fastq', 
								 'noTreat1':'CHe-XC-CON1_S5_L006_R1_001.fastq', 'noTreat2':'CHe-XC-CON2_S6_L006_R1_001.fastq', 
								 'kethoxal_rerun':'CHe-XC-MRE_S7_L006_R1_001.fastq', }

	library_info_dict['lib2']['fq'] = {'in_vivo_total_control_1':'CHe-XC-1_S1_L007_R1_001.fastq',
									   'in_vivo_total_control_2':'CHe-XC-2_S2_L007_R1_001.fastq',
									   'in_vivo_total_kethoxal_1':'CHe-XC-3_S3_L007_R1_001.fastq',
									   'in_vivo_total_kethoxal_2':'CHe-XC-4_S4_L007_R1_001.fastq',
									   'in_vivo_mRNA_control_1':'CHe-XC-5_S5_L007_R1_001.fastq',
									   'in_vivo_mRNA_control_2':'CHe-XC-6_S6_L007_R1_001.fastq',
									   'in_vivo_mRNA_kethoxal_1':'CHe-XC-7_S7_L007_R1_001.fastq',
									   'in_vivo_mRNA_kethoxal_2':'CHe-XC-8_S8_L007_R1_001.fastq',
									   'in_vitro_total_control_1':'CHe-XC-9_S9_L007_R1_001.fastq',
									   'in_vitro_total_control_2':'CHe-XC-10_S10_L007_R1_001.fastq',
									   'in_vitro_total_kethoxal_1':'CHe-XC-11_S11_L007_R1_001.fastq',
									   'in_vitro_total_kethoxal_2':'CHe-XC-12_S12_L007_R1_001.fastq',
									   'in_vitro_mRNA_control_1':'CHe-XC-13_S13_L007_R1_001.fastq',
									   'in_vitro_mRNA_control_2':'CHe-XC-14_S14_L007_R1_001.fastq',
									   'in_vitro_mRNA_kethoxal_1':'CHe-XC-15_S15_L007_R1_001.fastq',
									   'in_vitro_mRNA_kethoxal_2':'CHe-XC-16_S16_L007_R1_001.fastq',
									   }

	library_info_dict['lib3']['fq'] = {'in_vivo_mRNA_control':'CHe-Weng-C-1_S1_L006_R1_001.fastq',
									   'in_vivo_G4_control':'CHe-Weng-C-2_S2_L006_R1_001.fastq',
									   'in_vivo_Glucose_control':'CHe-Weng-C-3_S3_L006_R1_001.fastq',
									   'in_vivo_mRNA_kethoxal':'CHe-Weng-K-1_S4_L006_R1_001.fastq',
									   'in_vivo_G4_kethoxal':'CHe-Weng-K-2_S5_L006_R1_001.fastq',
									   'in_vivo_Glucose_kethoxal':'CHe-Weng-K-3_S6_L006_R1_001.fastq',
									   'in_vivo_t1min_kethoxal':'CHe-Weng-T-1_S7_L006_R1_001.fastq',
									   'in_vivo_t2.5min_kethoxal':'CHe-Weng-T-2_S8_L006_R1_001.fastq',
									   'in_vivo_t5min_kethoxal':'CHe-Weng-T-3_S9_L006_R1_001.fastq',
									   'in_vivo_t10min_kethoxal':'CHe-Weng-T-4_S10_L006_R1_001.fastq'}

	library_info_dict['lib6']['fq'] = {'noTreat_1':'C1_L4_A021.R1.fastq',
									   'noTreat_2':'C2_L4_A022.R1.fastq',
									   'kethoxal_1':'K1_L4_A024.R1.fastq',
									   'kethoxal_2':'K2_L4_A025.R1.fastq'}

	category_dict = {'rmdup':'rmdup.fastq', 'trimmed':'trimmed.fastq', 'sam':'sam', 'rpkm':'rpkm', 'rt':'rt', 'rRNA_sam':'rRNA.sam', 'rRNA_rt':'rRNA.rt'}

	for lib in ['lib5', 'lib2', 'lib3', 'lib6']:
		for category in ['rmdup', 'trimmed', 'sam', 'rpkm', 'rt', 'rRNA_sam', 'rRNA_rt']:
			for i,j in library_info_dict[lib]['fq'].items():
				library_info_dict[lib][category][i] = library_info_dict['dir_prefix']+'/'+library_info_dict[lib]['dir']+'/'+j.replace('fastq', category_dict[category])

	gj.print_dict(library_info_dict)

	return library_info_dict

def library_info():
	library_info_dict = nested_dict()
	library_info_dict['dir_prefix'] = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data'
	library_info_dict['lib']['dir'] = 'PE100'

	library_info_dict['lib']['fq'] = {'control1':'CHe-XC-M1C_S1_L005_R1_001.fastq', 'control2':'CHe-XC-M2C_S2_L005_R1_001.fastq', 
								 'kethoxal1':'CHe-XC-M1K_S3_L005_R1_001.fastq', 'kethoxal2':'CHe-XC-M2K_S4_L005_R1_001.fastq', 
								 'noTreat1':'CHe-XC-CON1_S5_L005_R1_001.fastq', 'noTreat2':'CHe-XC-CON2_S6_L005_R1_001.fastq', 
							         }
	library_info_dict['lib']['fq2'] = {'control1':'CHe-XC-M1C_S1_L005_R2_001.fastq', 'control2':'CHe-XC-M2C_S2_L005_R2_001.fastq', 
								 'kethoxal1':'CHe-XC-M1K_S3_L005_R2_001.fastq', 'kethoxal2':'CHe-XC-M2K_S4_L005_R2_001.fastq', 
								 'noTreat1':'CHe-XC-CON1_S5_L005_R2_001.fastq', 'noTreat2':'CHe-XC-CON2_S6_L005_R2_001.fastq', 
								 }

	category_dict = {'rmdup':'rmdup.fastq', 'trimmed':'trimmed.fastq', 'sam':'sam', 'rpkm':'rpkm', 'rt':'rt', 'rRNA_sam':'rRNA.sam'}

	for lib in ['lib']:
		for category in ['rmdup', 'trimmed', 'sam', 'rpkm', 'rt', 'rRNA_sam']:
			for i,j in library_info_dict[lib]['fq'].items():
				library_info_dict[lib][category][i] = library_info_dict['dir_prefix']+'/'+library_info_dict[lib]['dir']+'/'+j.replace('fastq', category_dict[category])

	gj.print_dict(library_info_dict)

	return library_info_dict

def read_len_dist(fq='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-M1K_S3_L005_R1_001.trimmed.fastq'):
	gj.printFuncRun('read_len_dist')
	gj.printFuncArgs()
	fq_len_txt = fq + '.len.txt'
	subprocess.call(["awk '{if(NR%%4==2) print length($1)}' %s| sort|uniq -c|sort -k2,2n > %s "%(fq, fq_len_txt)],shell=True) # use double % to escape 
	df = pd.read_csv(fq_len_txt, sep='\s+', header=None)
	df.columns = ['# of reads', 'read length']
	df_plot = df[['read length', '# of reads']]
	print df_plot
	fig, ax= plt.subplots(2,1, sharex=True)
	df.plot(ax=ax[0], x='read length', y='# of reads')
	df.plot(kind='scatter', ax=ax[0], x='read length', y='# of reads')

	df_trimlog = pd.read_csv(fq+'.trimlog', header=None, sep='\s+')
	df_trimlog.columns = ['seq_name', 'sample_name', 'survive_len', 'survive_start', 'survive_end', 'cut_len']
	df_trimlog = df_trimlog[df_trimlog['cut_len'] > 0]
	cut_len_ls = list(df_trimlog['cut_len'])

	n = [[i]*j for i,j in zip(df['read length'], df['# of reads'])]
	n = gj.ls_ls_flat(n)
	gj.cumulate_dist_plot([n, cut_len_ls],ls_ls_label=['kethoxal read length', 'kethoxal read cut length'],bins=40,title=None,ax=ax[1],savefn=None,xlabel='Length',ylabel=None,add_vline=None,add_hline=None,log2transform=0)

	plt.tight_layout()
	plt.savefig(fq+'.len.png')
	plt.close()

	gj.printFuncRun('read_len_dist')

def read_len_dist_all(savefn='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/read_len_dist_all.png'):
	gj.printFuncRun('read_len_dist_all')
	gj.printFuncArgs()
	library_info_dict = library_info()
	trimmed_dict = library_info_dict['lib']['trimmed']
	print trimmed_dict
	read_len_ls_ls = []
	read_cut_len_ls_ls = []
	fig, ax = plt.subplots(3,1, sharex=True, figsize=(14, 16))
	color_ls = gj.sns_color_ls()
	sample_ls = []
	for n, (i,j) in enumerate(trimmed_dict.items()):
		sample_ls.append(i)
		print i, j
		fq_len_txt = j + '.len.txt'
		trimlog = j + '.trimlog'
		df = pd.read_csv(fq_len_txt, sep='\s+', header=None)
		df.columns = ['# of reads', 'read length']
		df.plot(ax=ax[0], x='read length', y='# of reads', label=i)
		df_trimlog = pd.read_csv(trimlog, header=None, sep='\s+')
		df_trimlog.columns = ['seq_name', 'sample_name', 'survive_len', 'survive_start', 'survive_end', 'cut_len']
		df_trimlog = df_trimlog[df_trimlog['cut_len'] > 0]
		cut_len_ls = list(df_trimlog['cut_len'])
		n = [[i]*j for i,j in zip(df['read length'], df['# of reads'])]
		n = gj.ls_ls_flat(n)
		read_len_ls_ls.append(n)
		read_cut_len_ls_ls.append(cut_len_ls)
	gj.cumulate_dist_plot(read_len_ls_ls,ls_ls_label=sample_ls,bins=40,title=None,ax=ax[1],savefn=None,xlabel='Length',ylabel=None,add_vline=None,add_hline=None,log2transform=0)
	gj.cumulate_dist_plot(read_cut_len_ls_ls,ls_ls_label=sample_ls,bins=40,title=None,ax=ax[2],savefn=None,xlabel='Length',ylabel=None,add_vline=None,add_hline=None,log2transform=0)
	plt.tight_layout()
	plt.savefig(savefn)
	plt.close()
	gj.printFuncRun('read_len_dist_all')

def read_pair_len_dist(fastq1=None, fastq2=None, savefn=None):
	gj.printFuncRun('read_pair_len_dist')
	if fastq1 is None:
		fastq1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-M1K_S3_L005_R1_001.paired.fastqT'
	if fastq2 is None:
		fastq2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-M1K_S3_L005_R2_001.paired.fastqT'
	if savefn is None:
		savefn = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-M1K_S3_L005.paired.fastq.len.png'
	gj.printFuncArgs()
	read_len_ls1 = [] 
	with open(fastq1, 'r') as FQ1:
		for n,line in enumerate(FQ1):
			if n%4 == 1:
				read_len_ls1.append(len(line.strip()))
	read_len_ls2 = [] 
	with open(fastq2, 'r') as FQ2:
		for n,line in enumerate(FQ2):
			if n%4 == 1:
				read_len_ls2.append(len(line.strip()))
	df = pd.DataFrame({'read1':read_len_ls1, 'read2':read_len_ls2})
	print df.head()
	gj.df_sns_jointplot(col_str_x='read1',col_str_y='read2',savefn=savefn,df=df,list1='list1',list2='list2',xlim=None,ylim=None,x_y_lim_same=1,title_str='',title_suptitle='right',use_scale_x_y_lim=0,color=None,xlabel=None,ylabel=None)

	gj.printFuncRun('read_pair_len_dist')

def inner_distance(sam=None, output_prefix=None):
	gj.printFuncRun('inner_distance')
	gj.printFuncArgs()
	bed12 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10.transCoor.bed12'
	subprocess.call(["inner_distance.py -i %s -o %s -r %s"%(sam, output_prefix, bed12)],shell=True)
	gj.printFuncRun('inner_distance')

def FPKM_count(bam=None, output_prefix=None, rRNA=0):
	gj.printFuncRun('FPKM_count')
	gj.printFuncArgs()
	bed12 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10.transCoor.bed12'
	if rRNA:
		bed12 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10.rRNA.bed12'
	subprocess.call(["FPKM_count.py -i %s -o %s -r %s"%(bam, output_prefix, bed12)],shell=True)
	gj.printFuncRun('FPKM_count')

def main():
	"""
	parser = argparse.ArgumentParser(description='icshape processing pipeline')
	parser.add_argument('-fq', metavar='fq', help = 'input fastq file', type=str)
	parser.add_argument('-mapper', metavar='mapper', default='Bowtie2', choices=['bowtie2','bowtie','STAR'], help='mapping software', type=str)
	parser.add_argument('-index_dir', metavar='index_dir' ,help= 'reference fastq index', type=str)
	args = parser.parse_args()

	collapse_fq = args.fq.replace('fastq','rmdup.fastq')
	trimmed_fastq = args.fq.replace('fastq','trimmed.fastq')
	map_sam = args.fq.replace('fastq','sam')
	map_rRNA_sam = args.fq.replace('fastq','rRNA.sam')
	sam_rpkm = args.fq.replace('fastq','rpkm')
	sam_rt = args.fq.replace('fastq','rt')

	library_info_dict = library_info()
	"""
	library_info_dict = library_info2()
	#print library_info_dict

	#read_collapse(fastq=args.fq)
	#remove_adapter(fastq=collapse_fq)
	#mapping(fastq=trimmed_fastq)
	#rpkm_cal(sam=map_sam)
	#RT_cal(sam=map_sam)
	#rpkm_cal(sam=map_rRNA_sam)
	#RT_cal(sam=map_rRNA_sam)
	#library_info()

	#read_clean_map_rt(fastq=args.fq)
	#RT_correlation(rt1='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1C_S1_L006_R1_001.rt', rt2='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2C_S2_L006_R1_001.rt', rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control.corr.txt')
	#RT_correlation(rt1=library_info_dict['lib5']['rt']['kethoxal1'], rt2=library_info_dict['lib5']['rt']['kethoxal2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.corr.txt')
	#RT_correlation(rt1=library_info_dict['lib5']['rt']['kethoxal1'], rt2=library_info_dict['lib5']['rt']['kethoxal2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.corr.T2.b100.txt', coverage_cutoff=2, background_base_density=100)
	#RT_correlation(rt1=library_info_dict['lib5']['rt']['control1'], rt2=library_info_dict['lib5']['rt']['control2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control.corr.txt')
	#RT_correlation(rt1=library_info_dict['lib5']['rt']['control1'], rt2=library_info_dict['lib5']['rt']['control2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control.corr.T2.b100.txt', coverage_cutoff=2, background_base_density=100)
	#RT_correlation(rt1=library_info_dict['lib5']['rt']['noTreat1'], rt2=library_info_dict['lib5']['rt']['noTreat2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.corr.txt')
	#RT_correlation(rt1=library_info_dict['lib5']['rt']['noTreat1'], rt2=library_info_dict['lib5']['rt']['noTreat2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.corr.T2.b100.txt', coverage_cutoff=2, background_base_density=100)
	#RT_correlation(rt1=library_info_dict['lib5']['rt']['noTreat1'], rt2=library_info_dict['lib5']['rt']['noTreat2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.corr.T2.b100.txt', coverage_cutoff=2, background_base_density=100)
	#RT_combine(rt1=library_info_dict['lib5']['rt']['control1'], rt2=library_info_dict['lib5']['rt']['control2'], rt_comb='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control.rt')
	#RT_combine(rt1=library_info_dict['lib5']['rt']['kethoxal1'], rt2=library_info_dict['lib5']['rt']['kethoxal2'], rt_comb='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.rt')
	#RT_combine(rt1=library_info_dict['lib5']['rt']['noTreat1'],rt2=library_info_dict['lib5']['rt']['noTreat2'],rt_comb='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.rt')
	#RT_combine(rt1=library_info_dict['lib5']['rt']['control1'].replace('.rt','.rRNA.rt'), rt2=library_info_dict['lib5']['rt']['control2'].replace('.rt','.rRNA.rt'), rt_comb='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control_rRNA.rt')
	#RT_combine(rt1=library_info_dict['lib5']['rt']['kethoxal1'].replace('.rt','.rRNA.rt'), rt2=library_info_dict['lib5']['rt']['kethoxal2'].replace('.rt','.rRNA.rt'), rt_comb='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal_rRNA.rt')
	#RT_combine(rt1=library_info_dict['lib5']['rt']['noTreat1'].replace('.rt','.rRNA.rt'),rt2=library_info_dict['lib5']['rt']['noTreat2'].replace('.rt','.rRNA.rt'),rt_comb='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat_rRNA.rt')

	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control.rt')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.rt')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.rt')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control_rRNA.rt')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal_rRNA.rt')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat_rRNA.rt')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_control.rRNA.bam.RT')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_notreat.rRNA.bam.RT')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/control_control.rRNA.bam.RT')
	#RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/notreat_notreat.rRNA.bam.RT')

	#calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.normalized.rt', b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control.normalized.rt', icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq.tmp.out')
	#calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.normalized.rt', b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.normalized.rt', icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.tmp.out')
	#filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq.tmp.out')
	#filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.tmp.out')
	#calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100')

	#calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal_rRNA.normalized.rt', b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control_rRNA.normalized.rt', icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_rRNA.tmp.out')
	#calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal_rRNA.normalized.rt', b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat_rRNA.normalized.rt', icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_rRNA_noTreat.tmp.out')
	#filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_rRNA.tmp.out')
	#filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_rRNA_noTreat.tmp.out')

	# calc_enrich('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_control.rRNA.bam.normalized.RT', '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/control_control.rRNA.bam.normalized.RT', '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_control.icshape.tmp.out')
	# filter_enrich('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_control.icshape.tmp.out')
	# calc_enrich('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_control.rRNA.bam.normalized.RT', '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/notreat_notreat.rRNA.bam.normalized.RT', '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_notreat.icshape.tmp.out')
	# filter_enrich('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_notreat.icshape.tmp.out')
	# calc_enrich('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_control.rRNA.bam.normalized.RT2', '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/control_control.rRNA.bam.normalized.RT', '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_control.fulllength.tmp.out')
	# filter_enrich('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_control.fulllength.tmp.out')
	# calc_enrich('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_notreat.rRNA.bam.normalized.RT2', '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/notreat_notreat.rRNA.bam.normalized.RT', '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_notreat.fulllength.tmp.out')
	# filter_enrich('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/kethoxal_notreat.fulllength.tmp.out')

	"""
	x_test_ls = [i/float(100) for i in range(1,51)]
	for x in x_test_ls:
		calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal_rRNA.normalized.rt', b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control_rRNA.normalized.rt', icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/x_test/kethoxalseq_rRNA.tmp.out.'+str(x), x=x)
		calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal_rRNA.normalized.rt', b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat_rRNA.normalized.rt', icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/x_test/kethoxalseq_rRNA_noTreat.tmp.out.'+str(x), x=x)
		filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/x_test/kethoxalseq_rRNA.tmp.out.'+str(x))
		filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/x_test/kethoxalseq_rRNA_noTreat.tmp.out.'+str(x))
	"""

	fastq1 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-CON1_S5_L005_R1_001.fastq'
	fastq2 = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-CON1_S5_L005_R2_001.fastq'
	#read_collapse(fastq=fastq1)
	#remove_adapter_PE_new(fastq1=fastq1, fastq2=fastq2)
	#read_pair_len_dist()
	#remove_adapter(fastq='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-M2K_S4_L005_R1_001.fastq', trimmed_fastq='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-M2K_S4_L005_R1_001.trimmed.fastq')
	#read_len_dist(fq='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-M2K_S4_L005_R1_001.trimmed.fastq')
	#read_len_dist_all()
	#read_collapse_PE(fastq1, fastq2)
	#remove_adapter_PE(fastq1=fastq1.replace('fastq','rmdup.fastq'), fastq2=fastq2.replace('fastq','rmdup.fastq'))
	#mapping_PE(fastq1=fastq1.replace('fastq','trimmed.paired.fastq'), fastq2=fastq2.replace('fastq','trimmed.paired.fastq'))
	#mapping_PE(fastq1=fastq1.replace('fastq','paired.clip.fastq'), fastq2=fastq2.replace('fastq','paired.clip.fastq'))
	#gj.sam2sortedbam('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-M1C_S1_L005_R1_001.sam')
	#rpkm_cal(sam=library_info_dict['lib']['sam']['noTreat2'])
	#inner_distance(sam='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-CON2_S6_L005_R1_001.sam', output_prefix='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CON2')
	#FPKM_count(bam='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-CON_S6_L005_R1_001.rRNA.sorted.bam', output_prefix='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/PE100/CHe-XC-CON_S6_L005_R1_001.rRNA.sorted.bam', rRNA=1)

	##################################################################################
	# *** 16-11-14_7_library_total_Kethoxal_remove lib
	condition_ls = ['kethoxal1', 'kethoxal2', 'control1', 'control2', 'noTreat1', 'noTreat2']
	# for condition in condition_ls:
		# print "sam: %s"%(library_info_dict['lib5']['sam'][condition])
		# rpkm_cal(sam=library_info_dict['lib5']['sam'][condition])
		# RT_cal(sam=library_info_dict['lib5']['sam'][condition])
		# rpkm_cal(sam=library_info_dict['lib5']['sam'][condition].replace('.sam', '.map.plus.sam'))

		# print "rt: %s"%(library_info_dict['lib5']['rt'][condition])
		# rt = library_info_dict['lib5']['rt'][condition]
		# RT_normalize(rt=rt)
	condition_ls = ['kethoxal', 'control', 'noTreat']
	# for condition in condition_ls:
		# RT_combine(rt1=library_info_dict['lib5']['rt']['%s1'%(condition)], rt2=library_info_dict['lib5']['rt']['%s2'%(condition)], rt_comb='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/%s.rt'%(condition))
		# RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/%s.rt'%(condition))
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.normalized.rt', b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/control.normalized.rt', icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq.tmp.out')
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.normalized.rt', b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.normalized.rt', icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq.tmp.out')
	# # filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.tmp.out', average_coverage=1, background_base_density=20, skip_leading=5, skilp_tailing=30)
	# RT_combine(rt1='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.rfam.rt', 
	# 		   rt2='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON2_S6_L006_R1_001.rfam.rt', 
	# 		   rt_comb='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.rfam.rt')
	# RT_combine(rt1='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.rfam.rt', 
	# 		   rt2='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.rfam.rt', 
	# 		   rt_comb='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.rfam.rt')
	# RT_normalize(rt='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.rfam.rt')
	# RT_normalize(rt='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.rfam.rt')
	# calc_enrich(f_normalized_rt='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxal.rfam.normalized.rt',
	# 			b_normalized_rt='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/noTreat.rfam.normalized.rt',
	# 			icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.rfam.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/kethoxalseq_noTreat.rfam.tmp.out',
				  # average_coverage=1, background_base_density=10, skip_leading=5, skilp_tailing=30)
	##################################################################################

	##################################################################################
	# *** 16-08-08_16_library_invivo_invitro lib
	# condition_ls = ['in_vivo_total_control', 'in_vivo_total_kethoxal', 'in_vivo_mRNA_control', 'in_vivo_mRNA_kethoxal']
	# condition_ls = ['in_vitro_total_control', 'in_vitro_total_kethoxal', 'in_vitro_mRNA_control', 'in_vitro_mRNA_kethoxal']
	# for condition in condition_ls:
		# RT_correlation(rt1=library_info_dict['lib2']['rt'][condition+'_1'], rt2=library_info_dict['lib2']['rt'][condition+'_2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/%s.corr.txt'%(condition))
		# RT_correlation(rt1=library_info_dict['lib2']['rt'][condition+'_1'], rt2=library_info_dict['lib2']['rt'][condition+'_2'], rt_corr='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/%s.corr.T2.b100.txt'%(condition), coverage_cutoff=2, background_base_density=100)
		# RT_combine(rt1=library_info_dict['lib2']['rt'][condition+'_1'], rt2=library_info_dict['lib2']['rt'][condition+'_2'], rt_comb='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/%s.rt'%(condition))
		# RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/%s.rt'%(condition))
	"""
	condition_ls = ['in_vivo_total', 'in_vivo_mRNA', 'in_vitro_total', 'in_vitro_mRNA']
	for condition in condition_ls[0:]:
		# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/%s_kethoxal.normalized.rt'%(condition), 
		        # b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/%s_control.normalized.rt'%(condition), 
		        # icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/%s_kethoxal.tmp.out'%(condition))
		filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/%s_kethoxal.tmp.out'%(condition),average_coverage=2, background_base_density=20, skip_leading=5, skilp_tailing=30)
	"""
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.normalized.rt'%('in_vivo_G4'),
	# 	b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.normalized.rt'%('in_vivo_mRNA'), 
	# 	icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4OvermRNA_kethoxal.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4OvermRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4OvermRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=200, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4OvermRNA_kethoxal.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_G4OvermRNA_kethoxal.tmp.out',average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30)

	
	##################################################################################

	##################################################################################
	# *** 16-09-09_10_library_G4-glucose_time lib
	# condition_ls = ['in_vivo_mRNA', 'in_vivo_G4', 'in_vivo_Glucose']
	# for condition in condition_ls:
		# RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_control.rt'%(condition))
		# RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.rt'%(condition))
		# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.normalized.rt'%(condition), 
		#         b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_control.normalized.rt'%(condition), 
		#         icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.tmp.out'%(condition))
		# pass

	"""
	condition_ls = ['in_vivo_t1min', 'in_vivo_t2.5min', 'in_vivo_t5min', 'in_vivo_t10min']
	for condition in ['in_vivo_t2.5min']:
		print
		# RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.rt'%(condition))
		calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.normalized.rt'%(condition), 
		        b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-08-08_16_library_invivo_invitro/in_vivo_total_control.normalized.rt',  # all use in vivo total  control as background
		        icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.tmp.out'%(condition))
	"""
	# print library_info_dict['lib3']['rRNA_rt']
	# for sample,rRNA_sam in library_info_dict['lib3']['rRNA_sam'].items():
	# 	print "process: %s"%(rRNA_sam)
	# 	rpkm_cal(sam=rRNA_sam)
	# 	RT_cal(sam=rRNA_sam)
	# for sample,rRNA_rt in library_info_dict['lib3']['rRNA_rt'].items():
	# 	print "process: %s"%(rRNA_rt)
	# 	RT_normalize(rt=rRNA_rt)
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-K-1_S4_L006_R1_001.rRNA.normalized.rt',
	# 			b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-C-1_S1_L006_R1_001.rRNA.normalized.rt',
	# 			icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.rRNA.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.rRNA.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)

	# rpkm_cal(sam='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-K-1_S4_L006_R1_001.newrRNA.sam')
	# RT_cal(sam='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-K-1_S4_L006_R1_001.newrRNA.sam')
	# RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-K-1_S4_L006_R1_001.newrRNA.rt')
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-K-1_S4_L006_R1_001.newrRNA.normalized.rt',
	# 	b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/CHe-Weng-C-1_S1_L006_R1_001.newrRNA.normalized.rt',
	# 	icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.newrRNA.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/in_vivo_mRNA_kethoxal.newrRNA.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# condition_ls = ['in_vivo_mRNA', 'in_vivo_G4', 'in_vivo_Glucose'] #+ ['in_vivo_t1min', 'in_vivo_t2.5min', 'in_vivo_t5min', 'in_vivo_t10min']
	# for condition in condition_ls:
	# 	for T in [0,1,2]:
	# 		for t in [0,20,200]:
	# 			filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/%s_kethoxal.tmp.out'%(condition),average_coverage=T, background_base_density=t, skip_leading=5, skilp_tailing=30)
	##################################################################################

	"""
	library_info_dict['lib3']['fq'] = {'in_vivo_mRNA_control':'CHe-Weng-C-1_S1_L006_R1_001.fastq',
									   'in_vivo_G4_control':'CHe-Weng-C-2_S2_L006_R1_001.fastq',
									   'in_vivo_Glucose_control':'CHe-Weng-C-3_S3_L006_R1_001.fastq',
									   'in_vivo_mRNA_kethoxal':'CHe-Weng-K-1_S4_L006_R1_001.fastq',
									   'in_vivo_G4_kethoxal':'CHe-Weng-K-2_S5_L006_R1_001.fastq',
									   'in_vivo_Glucose_kethoxal':'CHe-Weng-K-3_S6_L006_R1_001.fastq',
									   'in_vivo_t1min_kethoxal':'CHe-Weng-T-1_S7_L006_R1_001.fastq',
									   'in_vivo_t2.5min_kethoxal':'CHe-Weng-T-2_S8_L006_R1_001.fastq',
									   'in_vivo_t5min_kethoxal':'CHe-Weng-T-3_S9_L006_R1_001.fastq',
									   'in_vivo_t10min_kethoxal':'CHe-Weng-T-4_S10_L006_R1_001.fastq'}
	"""

	##################################################################################
	# *** 19-04-03 lib, hela
	# condition_ls = ['in_vivo_mRNA', 'in_vivo_G4']
	# for condition in condition_ls:
	# 	RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/%s_control.rt'%(condition))
	# 	RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/%s_kethoxal.rt'%(condition))
	# 	calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/%s_kethoxal.normalized.rt'%(condition), 
	# 	        b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/%s_control.normalized.rt'%(condition), 
	# 	        icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/%s_kethoxal.tmp.out'%(condition))
	# for i in [1,2]:
	# 	for j in [100, 200]:
	# 		filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_mRNA_kethoxal.tmp.out',average_coverage=i, background_base_density=j, skip_leading=5, skilp_tailing=30)
	# 		filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4_kethoxal.tmp.out',average_coverage=i, background_base_density=j, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_mRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/%s_kethoxal.normalized.rt'%('in_vivo_G4'),
	# 	b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/%s_kethoxal.normalized.rt'%('in_vivo_mRNA'), 
	# 	icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4OvermRNA_kethoxal.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4OvermRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4OvermRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=200, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4OvermRNA_kethoxal.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-03/in_vivo_G4OvermRNA_kethoxal.tmp.out',average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30)

	# RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_kethoxal.rfam.rt')
	# RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_control.rfam.rt')
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_kethoxal.rfam.normalized.rt',
				# b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_control.rfam.normalized.rt',
				# icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_kethoxal.rfam.tmp.out')
	filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_kethoxal.rfam.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_kethoxal.rfam.tmp.out',average_coverage=1, background_base_density=200, skip_leading=5, skilp_tailing=30)
	filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_kethoxal.rfam.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-09-09_10_library_G4-glucose_time/Rfam/in_vivo_mRNA_kethoxal.rfam.tmp.out',average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30)
	##################################################################################

	##################################################################################
	# *** 19-04-19 lib, hela
	# condition_ls = ['in_vitro_mRNA', 'in_vitro_G4']
	# for condition in condition_ls:
	# 	RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/%s_control.rt'%(condition))
	# 	RT_normalize(rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/%s_kethoxal.rt'%(condition))
	# 	calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/%s_kethoxal.normalized.rt'%(condition), 
	# 	        b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/%s_control.normalized.rt'%(condition), 
	# 	        icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/%s_kethoxal.tmp.out'%(condition))
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vivo_mRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vivo_G4_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)

	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/%s_kethoxal.normalized.rt'%('in_vitro_G4'),
	# 	b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/%s_kethoxal.normalized.rt'%('in_vitro_mRNA'), 
	# 	icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4OvermRNA_kethoxal.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4OvermRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4OvermRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=200, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4OvermRNA_kethoxal.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4OvermRNA_kethoxal.tmp.out',average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30)

	# in vitro, rG4seq, KPDS as control
	# for i in ['repall']: #, 'rep1', 'rep2', 'rep3', 'rep4']:
		# RT_normalize(rt='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Kit_2016_NatureMethod_RG4/raw/merge/polyA-KPDS.%s.rt'%(i))
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4_kethoxal.normalized.rt', 
	# 	        b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4_controlrG4seqKPDS.normalized.rt', 
	# 	        icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4_kethoxal_controlrG4seqKPDS.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4_kethoxal_controlrG4seqKPDS.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4_kethoxal_controlrG4seqKPDS.tmp.out',average_coverage=1, background_base_density=200, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4_kethoxal_controlrG4seqKPDS.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_G4_kethoxal_controlrG4seqKPDS.tmp.out',average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30)

	# for i in ['repall']:
		# RT_normalize(rt='/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Kit_2016_NatureMethod_RG4/raw/merge_polyA-K/polyA-K.%s.rt'%(i))
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_mRNA_kethoxal.normalized.rt', 
	# 	        b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_mRNA_controlrG4seqK.normalized.rt', 
	# 	        icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_mRNA_kethoxal_controlrG4seqK.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_mRNA_kethoxal_controlrG4seqK.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_mRNA_kethoxal_controlrG4seqK.tmp.out',average_coverage=1, background_base_density=200, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_mRNA_kethoxal_controlrG4seqK.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-04-19/in_vitro_mRNA_kethoxal_controlrG4seqK.tmp.out',average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30)

	# 19-07-05, in vitro control(no treat)
	# for i in ['N1', 'N2', 'P1', 'P2']:
		# RT_normalize('/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/%s_combined_R2.rt'%(i))
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_kethoxal.normalized.rt', 
	# 	        b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_control.normalized.rt', 
	# 	        icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_kethoxal.tmp.out')
	# calc_enrich(f_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_kethoxal.normalized.rt', 
	# 	        b_normalized_rt='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_control.normalized.rt', 
	# 	        icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_kethoxal.tmp.out')
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_kethoxal.tmp.out',average_coverage=1, background_base_density=200, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_kethoxal.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_G4_kethoxal.tmp.out',average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_kethoxal.tmp.out',average_coverage=1, background_base_density=200, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_kethoxal.tmp.out',average_coverage=2, background_base_density=100, skip_leading=5, skilp_tailing=30)
	# filter_enrich(icshape_tmp_out='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/19-07-05/in_vitro_mRNA_kethoxal.tmp.out',average_coverage=2, background_base_density=200, skip_leading=5, skilp_tailing=30)
	##################################################################################

if __name__ == '__main__':
	main()
