import sys
import os 
import itertools 
import argparse
"""
2 input 
INPUT 1 - SAMPLE file - list of files that are finished downloading  
IMPUT 2 - TOTAL SAMPLE LIST - total list of samples that are in the project
OUTPUT - matching sample name 
"""
# EXTRACT metadat -> dict: key = SRR, value = rest of metadata 

def extract_done_list(sample_file):

	sample_list = []
	with open(sample_file) as samples:
		for sample in samples:
			sample_list.append(sample.rstrip())
	return sample_list

def extract_sample_name(total_sample_list):
	download_dict = {}
	count = 0
	# with open(download_file) as samples:
	# 	for sample in samples:
	# 		count = count + 1
	# 		download_dict[sample.split('/')[5].split('.bam')[0]] = sample
	# 		if count%100 == 0:
	# 			print sample.split('/')[5].split('.bam')[0]
	# return download_dict
	"""
	1. make dict for 
	"""
	with open(total_sample_list)


ap = argparse.ArgumentParser()
ap.add_argument('sample', help = 'list of files that are finished downloading')
ap.add_argument('download', help = 'The text file given by dbgap manifest')
ap.add_argument('output', help = 'the new metadata in csv format that includes the file name will be written on this output file')
args = ap.parse_args()

sample_list = extract_done_list(args.sample)
down_dict = extract_sample_name(args.download)

outfile = open(args.output, 'w')

count = 0
for key, value in down_dict.iteritems():
	if key not in sample_list:
		count += 1
		temp_line = value
		outfile.write(temp_line)
	else:
		print key

outfile.close()