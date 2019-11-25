import gj
import argparse
import subprocess
from nested_dict import nested_dict
import sys,os

import icshape

def get_dir_fastq(dir='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-03-25_library_of_concentration'):
	gj.printFuncRun('get_dir_fastq')
	gj.printFuncArgs()
	fn_ls = os.listdir(dir)
	fastq_fn_ls = [i for i in fn_ls if i.endswith('fastq') and 'Undetermined' not in i]
	print fastq_fn_ls
	fastq_fn_ls = [dir+'/'+i for i in fastq_fn_ls]
	gj.printFuncRun('get_dir_fastq')
	return fastq_fn_ls

def run_fastq(fq):
	run_fastq_log = fq+'.log'
	gj.printFuncRun('run_fastq')
	gj.printFuncArgs()

	LOG = open(run_fastq_log,'w')
	sys.stdout = LOG
	sys.stderr = LOG

	fq = fq.replace('trimmed.fastq','fastq')

	collapse_fq = fq.replace('fastq','rmdup.fastq')
	trimmed_fastq = fq.replace('fastq','trimmed.fastq')
	map_sam = fq.replace('fastq','sam')
	map_rRNA_sam = fq.replace('fastq','rRNA.sam')
	sam_rpkm = fq.replace('fastq','rpkm')
	sam_rt = fq.replace('fastq','rt')

	# for raw data
	icshape.read_collapse(fastq=fq)
	icshape.remove_adapter(fastq=collapse_fq)
	icshape.mapping(fastq=trimmed_fastq)
	icshape.rpkm_cal(sam=map_sam)
	icshape.RT_cal(sam=map_sam)
	icshape.rpkm_cal(sam=map_rRNA_sam)
	icshape.RT_cal(sam=map_rRNA_sam)

	# for clean data
	# icshape.mapping(fastq=trimmed_fastq)
	# icshape.rpkm_cal(sam=map_sam)
	# icshape.RT_cal(sam=map_sam)
	# icshape.rpkm_cal(sam=map_rRNA_sam)
	# icshape.RT_cal(sam=map_rRNA_sam)

	# for Rfam
	# map_rfam_sam = fq.replace('fastq','rfam.sam')
	# rfam_sam_rpkm = fq.replace('fastq','rfam.rpkm')
	# index = '/Share2/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/Rfam/Parsed_Structure/human.dot.dedup.fa'
	# icshape.map_rfam(fastq=trimmed_fastq, map_sam=map_rfam_sam, index=index)
	# icshape.rpkm_cal(sam=map_rfam_sam)
	# icshape.RT_cal(sam=map_rfam_sam)

	gj.printFuncRun('run_fastq')

	LOG.close()
	sys.stdout = sys.__stdout__

if __name__ == '__main__':
	#get_dir_fastq()
	run_fastq(fq=sys.argv[1])
