import sys
import os 
import itertools



### extract_sample_name
# IN = download file  each line has : gtex/exchange/GTEx_phs000424/sra_source_files/SRR1068687/G34712.GTEX-XXEK-0526-SM-4BRWD.1.bam.ncbi_enc 
# OUT = dict with - key = sample name including G code - G****.GTEX-[subject id]-[sample-id]......
def extract_sample_name(download_file):
	download_dict = {}
	with open(download_file) as samples:
		for sample in samples:
			download_dict[sample.split('/')[5].split('.bam')[0]] = sample
	return download_dict


### extract done files 
# IN = done sample file - each lien contains sample that is finished
# OUT = list of names - G****.GTEX-[subject ID] .... 
def extract_done(done_file):
	done_list = []
	with open(done_file) as samples:
		for sample in samples:
			sample_name = sample.split('.')[0:3]
			done_list.append(sample_name)
	return done_list	

ap = argparse.ArgumentParser()
ap.add_argument('download_manifest', help = 'Download manifest file - i.e. gtex/exchange/GTEx_phs000424/sra_source_files/SRR1068687/G34712.GTEX-XXEK-0526-SM-4BRWD.1.bam.ncbi_enc')
ap.add_argument('done_file', help = 'List of files that are done - i.e. one line is : G12223.GTEx-............')
ap.add_argument('out_File', help = 'OUTPUT will be generated in this file - it will have the things that needs to be re-downloaded and run')
args = ap.parse_args()

down_dict = extract_sample_name(ap.download_manifest)
done_sample_list = extract_done(ap.done_file)
with open(ap.out_File) as output: 
	for key, value in down_dict.iteritems():
		if key not in done_sample_list:
			output.write(str(value + '\n'))
		else:
			print "%s is found!" % (key)