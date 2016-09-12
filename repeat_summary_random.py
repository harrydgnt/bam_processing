import sys
import os 
import itertools
import argparse
import numpy as np 
import pandas as pd 
import time


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
class Dictlist(dict):
    def __setitem__(self, key, value):
        try:
            self[key]
        except KeyError:
            super(Dictlist, self).__setitem__(key, [])
        self[key].append(value)

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
			element_list.append(str(element.split('>')[1].rstrip().split()[0].split('/')[0]))
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

def extract_reads_temp(element_list, repeat_file):
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
		status = 1
		count = -1

		# status 0 = new read
		# status 1 = old read - same as the previous one 
		current_read_name = ""
		current_read_highest_score = -1
		current_read_highest_element = ""  
		previous_read_name = ""
		previous_score = -1 
		previous_element = "" 
		current_element_list = []
		for line in lines:
			try:
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
				element = str(line.split()[1])

				# ORIGINAL - no filtering just number
				# if name != current_read_name:
				# 	# ADD THE NUMBER 
				# 	element = str(line.split()[1])
				# 	current_read_name = name
				# 	element_dict[element] += 1
				# 	status = 0
				# 	num_reads += 1
				# 	current_read_highest_score = score
				# else: # if the name is the same
				# 	if status == 0 and score >= current_read_highest_score:
				# 		status = 1 
				# 		num_multimapped += 1
				# 	else:
				# 		continue


				# # ELIMINATE THE DUPLICATE ONES 
				# if name != previous_read_name:
				# 	# # ADD THE NUMBER 
				# 	# element = str(line.split()[1])
				# 	# current_read_name = name
				# 	# element_dict[element] += 1
				# 	# status = 0
				# 	# num_reads += 1
				# 	# current_read_highest_score = score
				# 	if status == 0:
				# 		element_dict[previous_element] += 1
					 
				# 	previous_read_name = name
				# 	previous_score = score
				# 	previous_element = str(line.split()[1])
				# 	status = 0
				# else: # if the name is the same
				# 	if status == 0 and score >= current_read_highest_score:
				# 		status = 1 
				# 		num_multimapped += 1
				# 	else:
				# 		continue

				# # RANDOMLY SELECT ONE OR THE OTHER 
				# if name != previous_read_name: 
				# 	num_reads += 1
				# 	if count != -1:
				# 		# select random one
				# 		selected_element = current_element_list[np.random.randint(1000) % len(current_element_list)]
				# 		# np.random.randint()
				# 		element_dict[selected_element] += 1 # add one for that element 
				# 		del current_element_list[:]
				# 	current_element_list.append(element)
				# 	current_read_highest_score = score
				# 	status = 0
				# 	count += 1
				# else:
				# 	if status == 0 and score >= current_read_highest_score:
				# 		current_element_list.append(element)
				# 		num_multimapped += 1
				# 	else:
				# 		continue

				# RANDOMLY DISTRIBUTE





			except IndexError:
				print "index error: ", line



	return element_dict, num_reads, num_multimapped


def extract_reads(element_list, repeat_file, ):
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
		status = 1
		count = -1

		# status 0 = new read
		# status 1 = old read - same as the previous one 
		current_read_name = ""
		current_read_highest_score = -1
		current_read_highest_element = ""  
		previous_read_name = ""
		previous_score = -1 
		previous_element = "" 
		current_element_list = []
		for line in lines:
			try:

				# below is assmuing two things: 
				# 1. we use the first entry 
				# 2. highest score element is on top 
				name = str(line.split()[0])
				score = float(line.split()[11])
				element = str(line.split()[1])

				if current_read_name != name: # different read 
					current_read_name = name
					current_read_highest_score = score
					if status != 1:
						for item in current_element_list:
							element_dict[item] += float(1/len(current_element_list))
					del current_element_list[:]
					current_element_list.append(element)
				else: # if same read
					if current_read_highest_score == score:
						current_element_list.append(element)
					elif current_read_highest_score > score:
						continue
					else: 
						print "ERROR - not the highest"

			except IndexError:
				print "index error: ", line



	return element_dict, num_reads, num_multimapped


def extract_bed(bed_file, element_list):
	"""
	bed 7 for repeat makser 
	[chr,start,end,strand, element, class, family]

	IN - repeat masker bed file 
	OUT = dict? list of things?
	elements have multiple positions across the genome 
	dict(dict(int)) - dict of (element -> position)
	"""
	# element_dict = {} # this is used for count
	element_dictlist = Dictlist() 
	for line in bed_file:
		try:
			element = line.split()[4]
			# element_dict[element] = 0
			if element in element_list:
				element_dictlist[element] = line
			else:
				continue
		except IndexError:
			print line
	return element_dictlist

def make_merge_dataframe(original_df, dict_to_add, sample_name):
	temp_df = pd.DataFrame.from_dict(dict_to_add, orient="index")
	temp_df.columns = [sample_name]
	return pd.merge(original_df, temp_df, left_index = True, right_index = True)

def edit_dict(input_dict, position_dict_list):
	new_dict = {}
	for element, num_reads in input_dict.iteritems:
		for item in position_dict_list[element]:
			new_dict[item] += float(num_reads/len(position_dict_list[element]))
	return new_dict

def main(element_list_file, sample_dir, sample_list, outfile, bed_file):
	start_time = time.time()
	element_list = extract_element(element_list_file)
	summary_df = pd.DataFrame()
	os.chdir(sample_dir)
	count = 0
	num_processed = 0
	pos_dict = extract_bed(bed_file,element_list)
	with open(sample_list, 'r') as samples:
		for sample in samples:
			num_processed += 1
			sample = sample.rstrip()
			temp_dict, num_reads, num_multimapped = extract_reads(element_list, sample)
			# edit step
			current_dict = edit_dict(temp_dict, pos_dict)

			sample_name = '.'.join(sample.split('.')[1:3])
			# print "Processing: ", sample_name
			print num_reads, num_multimapped
			if count == 0:
				summary_df = pd.DataFrame.from_dict(current_dict, orient = "index")
				summary_df.columns = [sample_name]
				count = -1
			else: 
				summary_df = make_merge_dataframe(summary_df, current_dict, sample_name)
			if num_processed % 1000 == 0:
				print "num_processed : ", num_processed
				print "elapsed_time = ", time.time() - start_time
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
test_sample = 'sample.txt'
outfile = '/u/home/s/serghei/result/repeat_summary_random.txt'
bedfile = '/u/home/h/harryyan/project-eeskin/gtex_repeat/db/library/repetitiveElements_hg19.bed7'
main(element_file,test_dir, test_sample, outfile, bedfile)


