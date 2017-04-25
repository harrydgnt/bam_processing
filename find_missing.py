import sys
import os 
import itertools
import argparse




##### 
# MATCH SRR 
###########

def import_done(input_file):
	"""
	input = files with line such as :
	C1422.GTEX-S341-0001.1.unmapped_after_rRNA_lostHuman.fasta.gz
	G16668.GTEX-PX3G-0126.3.unmapped_after_rRNA_lostHuman.fasta.gz
	First three entries split by . used for ID 
	OUTPUT - list of finished items
	"""
	finished_item_list = []
	count = 0
	with open(input_file, 'r') as f:
		for line in f:
			sample = ".".join(line.split('.')[0:3])
			finished_item_list.append(sample)
			
			count += 1
			if count % 1000 == 0: 
				print "finished importing %i completed samples" % (count)
	return finished_item_list

def get_SRR_dict(input_manifest_file):
	"""
	INPUT = files with lines:
	gtex/exchange/GTEx_phs000424/sra_source_files/SRR599346/G16763.GTEX-P78B-0526.4.bam.ncbi_enc
	gtex/exchange/GTEx_phs000424/sra_source_files/SRR598148/G16584.GTEX-NPJ8-0426.2.bam.ncbi_enc
	OUTPUT = dict : key -> sample name, value -> SRR
	"""
	SRR_dict = {}
	count = 0
	with open(input_manifest_file, 'r') as f:
		for line in f:
			sample_name = '.'.join(line.split('/')[-1].split('.')[0:3])
			# print sample_name
			SRR_number = line.split('/')[-2].split('_')[0]
			SRR_dict[sample_name] = SRR_number

			count+=1
			if count % 1000 == 0:
				print "finished importing manifest for %i samples" % (count)
	return SRR_dict


def get_existing_SRR(completed_samples, SRR_dict):
	completed_SRR_list = []
	for sample in completed_samples:
		try:
			completed_SRR_list.append(SRR_dict[sample])
		except KeyError:
			print "%s failed to find matching SRR" % (sample)
	return completed_SRR_list

def get_manifest_dict(input_manifest_file):
	manifest_dict = {}
	count = 0
	with open(input_manifest_file, 'r') as f:
		for line in f:
			sample_name = '.'.join(line.split('/')[-1].split('.')[0:3])
			# print sample_name
			SRR_number = line.split('/')[-2].split('_')[0]
			manifest_dict[SRR_number] = line.rstrip()
			count+=1
			if count % 1000 == 0:
				print "finished importing manifest for %i samples" % (count)
	return manifest_dict

def get_incomplete_manifest(completed_SRR_list, manifest_dict, output_manifest):
	count = 0
	count_incomplete = 0
	with open(output_manifest, 'w') as outfile:

		for key, value in manifest_dict.iteritems():
			count += 1
			if key not in completed_SRR_list:
				count_incomplete += 1
				outfile.write(value + '\n')
	print count, count_incomplete


def test():
	input_done = './done_preliminary_samples.txt' 
	input_manifest = './gtex_v6_bam_manifest.csv'
	output_manifest = './incomplete_manifest.txt'
	completed_sample_list = import_done(input_done)
	#print completed_sample_list
	SRR_dict = get_SRR_dict(input_manifest)
	# print SRR_dict
	completed_SRR_list = get_existing_SRR(completed_sample_list, SRR_dict)
	# print completed_SRR_list
	manifest_dict = get_manifest_dict(input_manifest)
	# print manifest_dict
	get_incomplete_manifest(completed_SRR_list, manifest_dict, output_manifest)
if __name__ == "__main__":
	print "x"
	test()