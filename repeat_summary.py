import sys
import os 
import itertools
import argparse
import numpy as np 
import pandas as pd 

####### 
# DOCUMENTATION
######

# list of elements available at  ~/project-eeskin/gtex_repeat/repeat_elements.txt
# 
# Example line from TSV output : http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# C37B2ACXX140318:7:2116:19986:26481/1    AluSg4-SINE1/7SL-Primates       97.37   76      2       0       1       76      205     280     4e-32    130
# read_name 	repeat_elements 	percent_identity	alignment_length	num_mismatch	num_gap_opening		qstart		qend	ref_start	ref_end		expect_value	bit_score
#
# Serghei's cutoff
# if eValue<1e-05 and alignmentLength>=0.8*readLength and identity>=0.9*readLength:
# 
#
#
# TODO - what to do with the multi mapped reads 
# 
# 
# 
# 
# 
# 
# 
# 
"""
PLAN
1. extract element names from the list of elements
- possibly dict(list) - key = element name, value = list 
--> IN = list of elements
--> OUT = element_dict

* RECURSIVE 
	2. extract each read's element
	- use the cutoff?
	- what to do with the multi-mapped read?  
	--> IN = element_dict, tsv file 
	--> OUT = element_dict with updated num_reads/element

3. Make matrix
- stack csv files? 
- everytime (2) is run?
- is it going to be too big?
	- if so write it on txt and convert it after? 
	- appending doesnt take much memory 

4. Transpose the matrix 

"""
def extract_element(element_list_file):
	'''
	1. extract element names from the list of elements
	- possibly dict(list) - key = element name, value = list 
	--> IN = list elements file
	--> OUT = element_dict
	'''	
	element_list = []
	with open(element_list_file, 'r') as elements:
		for element in elements:
			element_list.append(str(element.split('>')[1].rstrip().split()[0]))
	return element_list

# def make_element_dict(element_list):
# 	'''
# 	1. extract element names from the list of elements
# 	- possibly dict(list) - key = element name, value = list 
# 	--> IN = list of elements
# 	--> OUT = element_dict
# 	'''	
# 	element_dict = {}
# 	with open(element_list_file, 'r') as elements:
# 		for element in elements:
# 			element_dict[str(element.split('>')[1].rstrip().split()[0])] = 0
# 			# print str(element.split('>')[1].rstrip().split()[
# 	return element_dict	

def extract_reads(element_list, repeat_file):
	'''
	2. extract each read's element
	- use the cutoff?
	- what to do with the multi-mapped read?  
	--> IN = element_dict, tsv file 
	--> OUT = element_dict with updated num_reads/element
	'''

	# make dict from element list 
	element_dict = {}
	for item in element_list:
		element_dict[item] = 0
	num_multimapped = 0
	num_reads = 0
	# add number
	with open(repeat_file, 'r') as lines:
		status = 0;

		# status 0 = new read
		# status 1 = old read - same as the previous one 
		current_read_name = ""
		current_read_highest_score = -1
		current_read_highest_element = ""  
		for line in lines:
			# if status == 0:
			# 	current_read_name = str(line.split()[0])
			# 	current_read_highest_score = float(line.split()[11])
			# 	current_read_highest_element = str(line.split()[1])


			# 	status = 1
			# elif status == 1: 
			# 	name = str(line.split([0]))
			# 	score = int(line.split()[11])
			# 	if name == current_read_name: # if we are looking at the same read
			# 		continue
			# 	else: # if it is different read
			# 		status = 0




			# below is assmuing two things: 
			# 1. we use the first entry 
			# 2. highest score element is on top 
			name = str(line.split()[0])
			score = float(line.split()[11])

			if name != current_read_name:
				# ADD THE NUMBER 
				element = str(line.split()[1])
				current_read_name = name
				element_dict[element] += 1
				status = 0
				num_reads += 1
			else: # if the name is the same
				if status == 0:
					status = 1 
					num_multimapped += 1
				else:
					continue
	return element_dict, num_reads, num_multimapped

def make_merge_dataframe(original_df, dict_to_add, sample_name):
	temp_df = pd.DataFrame.from_dict(dict_to_add, orient="index")
	temp_df.columns = [sample_name]
	return pd.merge(original_df, temp_df, left_index = True, right_index = True)

def main(element_list_file, sample_dir, sample_list, outfile):
	element_list = extract_element(element_list_file)
	summary_df = pd.DataFrame()
	os.chdir(sample_dir)
	count = 0
	with open(sample_list, 'r') as samples:
		for sample in samples:
			sample = sample.rstrip()
			current_dict, num_reads, num_multimapped = extract_reads(element_list, sample)
			sample_name = '.'.join(sample.split('.')[1:3])
			# print "Processing: ", sample_name
			print num_reads, num_multimapped
			if count == 0:
				summary_df = pd.DataFrame.from_dict(current_dict, orient = "index")
				summary_df.columns = [sample_name]
				count = -1
			else: 
				summary_df = make_merge_dataframe(summary_df, current_dict, sample_name)
	summary_df.to_csv(outfile, sep = "\t")

	# TODO - after indexing everything, add class. family information

def test():
	element_file = '/u/home/h/harryyan/project-eeskin/gtex_repeat/repeat_elements.txt'
	test_file = '/u/home/h/harryyan/project-eeskin/gtex_repeat/G60826.GTEX-13O3O-0011-R1b.2.unmapped_after_rRNA_lostHuman.fasta_lostRepeats_blastFormat6.tsv'
	test_file_two = '/u/home/h/harryyan/project-eeskin/gtex_repeat/G20897.GTEX-SN8G-0006-SM-32PLD.1.unmapped_lostRepeats_blastFormat6.tsv'
	elem_list = extract_element(element_file)
	new_dict = extract_reads(elem_list, test_file)
	new_dict_two = extract_reads(elem_list, test_file_two)
	# for key,value in new_dict.iteritems():
	# 	print key, "\t", value

	test_df = pd.DataFrame()
	test_df = pd.DataFrame.from_dict(new_dict,orient="index")
	test_df.columns = ["G60826"]
	print test_df

	test_two = pd.DataFrame.from_dict(new_dict_two,orient="index")
	print test_two

	test_merge_df = pd.merge(test_df,test_two, left_index=True, right_index=True)
	print test_merge_df	

	print "shape, ", test_df.shape, test_two.shape, test_merge_df.shape
	# TODO - ADD sample name for each column 


# test()

element_file = '/u/home/h/harryyan/project-eeskin/gtex_repeat/repeat_elements.txt'
test_dir = '/u/home/s/serghei/result/unmapped/repeat_unmapped/'
test_sample = 'test_300.sample'
outfile = '/u/home/s/serghei/result/repeat_summary.txt'
main(element_file,test_dir, test_sample, outfile)


