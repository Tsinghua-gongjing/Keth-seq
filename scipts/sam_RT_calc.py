import subprocess
from pyfasta import Fasta
import gj
import sys

def bamToBed(bam=None):
	if bam is None:
		bam = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.bam'
	bed = bam.replace('.bam', '.bed')
	subprocess.call(["bedtools bamtobed -i %s > %s"%(bam, bed)], shell=True)

	return bed

def samToBam(sam=None):
	if sam is None:
		sam = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M1K_S3_L006_R1_001.sam'
	bam = sam.replace('.sam', '.bam')
	subprocess.call(["samtools view -b -S %s > %s"%(sam, bam)], shell=True)

	return bam

def sortBamIndex(bam):
	# sort_bam_prefix = bam.replace('.bam', '.sorted.bam')
	sort_bam = bam.replace('.bam', '.sorted.bam')
	subprocess.call(["samtools sort %s -o %s; samtools index %s"%(bam, sort_bam, sort_bam)], shell=True)

def extract_sam_map_plus(sam=None):
	if sam is None:
		sam = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-M2K_S4_L006_R1_001.sam'
	map_sam = sam.replace('.sam', '.map.sam')
	map_plus_sam = map_sam.replace('.sam', '.plus.sam')
	subprocess.call(["grep ^@ %s > %s ; samtools view -F 0x04 %s >> %s"%(sam, map_sam, sam, map_sam)], shell=True)
	subprocess.call(["grep ^@ %s > %s ; samtools view -F 0x10 %s >> %s"%(map_sam, map_plus_sam, map_sam, map_plus_sam)], shell=True)
	map_plus_bam = samToBam(sam=map_plus_sam)
	sortBamIndex(map_plus_bam)
	map_plus_bed = bamToBed(bam=map_plus_bam)


def read_fa(fa='/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/mm10/transcriptome/mm10_transcriptome.fa'):
	gj.printFuncRun('read_fa')
	gj.printFuncArgs()
	fa_dict1 = Fasta(fa, key_fn=lambda key:key.split("\t")[0])
	fa_dict = {i.split()[0].split('.')[0]:j[0:] for i,j in fa_dict1.items()}
	print fa_dict.keys()[0:3]
	gj.printFuncRun('read_fa')
	return fa_dict

def ol_bed_filter(ol_bed=None):
	if ol_bed is None:
		ol_bed = '/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/M1K_S3-M1C_S1.ol.bed'
	savefn = ol_bed.replace('.bed', '.filter.bed')
	SAVEFN = open(savefn, 'w')
	fa_dict = read_fa()
	with open(ol_bed, 'r') as BED:
		for n,line in enumerate(BED):
			if n % 1000000 == 0:
				print "process: %s"%(n)
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			# arr = line.split('\t')
			tx1,tx1_start,tx1_end,read1,score1,tx1_strand,tx2,tx2_start,tx2_end,read2,score2,tx2_strand,ol_base_num = line.split('\t')
			if tx1_end != tx2_end:
				continue
			if tx1_start == tx2_start:
				continue
			try:
				RT_base = fa_dict[tx1.split('.')[0]][int(tx1_end)]
			except:
				continue
			if RT_base != 'G':
				continue
			print >>SAVEFN, line
	SAVEFN.close()

def main():
	# samToBam()
	# bamToBed()
	# read_fa()
	# ol_bed_filter()
	# extract_sam_map_plus()
	sortBamIndex('/Share/home/zhangqf5/gongjing/Kethoxal_RNA_structure/data/16-11-14_7_library_total_Kethoxal_remove/CHe-XC-CON1_S5_L006_R1_001.map.plus.sorted.bam')

if __name__ == '__main__':
	# main()
	extract_sam_map_plus(sam=sys.argv[1])

"""
grep ^@ CHe-XC-M1C_S1_L006_R1_001.sam > CHe-XC-M1C_S1_L006_R1_001.map.sam # header
samtools view -F 0x04 CHe-XC-M1C_S1_L006_R1_001.sam >> CHe-XC-M1C_S1_L006_R1_001.map.sam # output mapped read

sam2bed CHe-XC-M1C_S1_L006_R1_001.map.sam CHe-XC-M1C_S1_L006_R1_001.map.bed
grep "+" CHe-XC-M1C_S1_L006_R1_001.map.bed > CHe-XC-M1C_S1_L006_R1_001.map.plus.bed # keep plus strand

bedtools intersect -a CHe-XC-M1K_S3_L006_R1_001.map.plus.bed -b CHe-XC-M1C_S1_L006_R1_001.map.plus.bed -wo -f 1 -s > M1K_S3-M1C_S1.ol.bed
"""