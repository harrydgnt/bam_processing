import sys
import os 
import itertools 
import argparse
"""
2 input 
INPUT 1 - SAMPLE file - list of files that are finished downloading  
IMPUT 2 - Download location manifest
	i.e. gtex/exchange/GTEx_phs000424/sra_source_files/SRR1068687/G34712.GTEX-XXEK-0526-SM-4BRWD.1.bam.ncbi_enc
OUTPUT - matching sample name 
"""
# EXTRACT metadat -> dict: key = SRR, value = rest of metadata 

def extract_done_list(sample_file):

	sample_list = []
	with open(sample_file) as samples:
		for sample in samples:
			sample_list.append(sample)
	return sample_list

def extract_sample_name(download_file):
	download_dict = {}
	with open(download_file) as samples:
		for sample in samples:
			download_dict[sample.split('/')[5].split('.bam')[0]] = sample
	return download_dict



ap = argparse.ArgumentParser()
ap.add_argument('sample', help = 'list of files that are finished downloading')
ap.add_argument('download', help = 'The text file given by dbgap manifest')
ap.add_argument('output', help = 'the new metadata in csv format that includes the file name will be written on this output file')
args = ap.parse_args()

sample_list = extract_done_list(args.sample)
down_dict = extract_sample_name(args.download)

outfile = open(args.output, 'w')

for key, value in down_dict.iteritems():
	if key not in sample_list:
		temp_line = value
		outfile.write(temp)

outfile.close()