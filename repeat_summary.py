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
	--> IN = list of elements
	--> OUT = element_dict
	'''	
	element_dict = {}
	with open(element_list_file, 'w') as elements:
		for element in elements:
			element_dict[element] = 0
	return element_dict	

def extract_reads(element_dict, repeat_file):
	'''
	2. extract each read's element
	- use the cutoff?
	- what to do with the multi-mapped read?  
	--> IN = element_dict, tsv file 
	--> OUT = element_dict with updated num_reads/element
	'''
	with open(repeat_file, 'w') as lines:
		status = 0 
		# status 0 = new read
		# status 1 = old read - same as the previous one 
		current_read_name = ""
		current_read_highest_score = -1
		current_read_highest_element = ""  
		for line in lines:
			# if status == 0:
			# 	current_read_name = str(line.split()[0])
			# 	current_read_highest_score = int(line.split()[11])
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
			name = str(line.split([0]))
			score = int(line.split()[11])

			if name != current_read_name:
				# ADD THE NUMBER 
				element = str(line.split()[1])
				current_read_name = name
				element_dict[element] += 1
			else: # if the name is the same
				continue
	return element_dict	


def test():
	element_file = '/u/home/h/harryyan/project-eeskin/gtex_repeat/db/repeat_elements.txt'
	test_file = '/u/home/h/harryyan/project-eeskin/gtex_repeat/G60826.GTEX-13O3O-0011-R1b.2.unmapped_after_rRNA_lostHuman.fasta_lostRepeats_blastFormat6.tsv'
	elem_dict = extract_element(element_file)
	new_dict = extract_reads(elem_dict, test_file)
	for key,value in new_dict.iteritems():
		print key, "\t", value