import sys
import os 
import itertools 
import argparse
"""
2 input 
INPUT 1 - SAMPLE Manifest - SRR number on the first column and the rest in the back in CSV format 
IMPUT 2 - Download location manifest
	i.e. gtex/exchange/GTEx_phs000424/sra_source_files/SRR1068687/G34712.GTEX-XXEK-0526-SM-4BRWD.1.bam.ncbi_enc
OUTPUT - matching sample name 
"""


# EXTRACT metadat -> dict: key = SRR, value = rest of metadata 

def extract_metadata(sample_file):

	sample_dict = {}
	with open(sample_file) as samples:
		for sample in samples:
			sample_dict[sample.split(',')[0]] = sample.split(',')[1:]
	return sample_dict

def extract_sample_name(download_file):
	download_dict = {}
	with open(download_file) as samples:
		for sample in samples:
			download_dict[sample.split('/')[4]] = sample.split('/')[5].split('.bam')[0]

ap = argparse.ArgumentParser()
ap.add_argument('metadata', help = 'SRR####### on the first column in CSV format')
ap.add_argument('download', help = 'The text file given by dbgap manifest')
args = ap.parse_args()

print args