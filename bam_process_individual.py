import sys
import os
import argparse 
import numpy as np 



###############################################################
# Arguments
###############################################################
ap = argparse.ArgumentParser()
necessary_arguments = ap.add_argument_group('Necessary Inputs')
necessary_arguments.add_argument('ncbi_enc_dir', help="dir with *.ncbi_enc file")
necessary_arguments.add_argument('ncbi_enc_file', help="ncbi_enc file")
necessary_arguments.add_argument('decryption_depository', help="dir with decryption key installed, including vdb-decrypt")
necessary_arguments.add_argument('dump_dir', help="dir where decrypted bam")
necessary_arguments.add_argument('gene_list_file', help="list file of the genes to be extracted in bam format")
necessary_arguments.add_argument('immune_gene_list_file', help="list file of the immune genes to be extracted in bam format")

analysis_arguments = ap.add_argument_group('Analysis Options')
analysis_arguments.add_argument('--htseq', help="HTSeq-Count for gene count")
analysis_arguments.add_argument('--unmap_count', help = "Extract num_unmapped and num_total for the sample")
analysis_arguments.add_argument('--repeat', help = "Repeat Profile for mapped reads")
analysis_arguments.add_argument('--gene', help = "Extract genes")
analysis_arguments.add_argument('--immune', help = "Extract immune genes in bam format")
analysis_arguments.add_argument('--MT', help = "Extract MT reads")

optional_arguments = ap.add_argument_group('Optional Arguments')
optional_arguments.add_argument('--dev', help = "Dev option - keep the unmapped bam file")
optional_arguments.add_argument('--qsub', help = "The command will be submitted through qsub automatically")
#maybe running directly in python? 

args = ap.parse_args()


###############################################################
# Analysis Argument Handling 
###############################################################

if not args.htseq and not args.unmap_count and not args.repeat and not args.gene and not args.immune and not args.MT:
	args.htseq = True
	args.unmap_count = True
	args.repeat = True
	args.gene = True
	args.immune = True
	args.MT = True



###############################################################
# Source Code 
###############################################################

# I/O operation

try:
	os.mkdir(args.dump_dir)

os.chdir(args.dump_dir)

item = str(args.ncbi_enc_file.split('.ncbi_enc')[0])
item_name = str(item.split('.')[0])
ext_item_name = str(item.split('.bam')[0])

os.mkdir(ext_item_name)
os.chdir(ext_item_name)

command_output_file = "process_" + ext_item_name

#with open(command_output_file, 'a') as cmd_out:

# Modules
cmd = "#!/bin/bash \n"
cmd += ". /u/local/Modules/default/init/modules.sh \n"
cmd += "module load samtools \n"
cmd += "module load bamtools \n"

# Decryption 

cmd += "mv %s/%s %s \n" % (args.ncbi_enc_dir, args.ncbi_enc_file, args.decryption_depository)
cmd += "cd %s \n" % (args.decryption_depository)
cmd += "%s/vdb-decrypt -f %s/%s \n" % (args.decryption_depository, args.decryption_depository, args.ncbi_enc_file)
cmd += "mv %s/%s %s/%s \n" % (args.decryption_depository, item, args.dump_dir, ext_item_name)
cmd += "cd %s/%s \n" % (args.dump_dir, ext_item_name)

# Indexing the bam

cmd += "samtools index %s/%s/%s \n" % (args.dump_dir, ext_item_name, item)

# HTSeq count

if args.htseq:
	gtf = "/u/home/s/serghei/project/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
	cmd += "samtools sort -n %s/%s/%s %s/%s/%s \n" % (args.dump_dir, ext_item_name, item, args.dump_dir, ext_item_name, item)
	cmd += "python /u/home/h/harryyan/project-eeskin/utilities/HTSeq-0.6.1/scripts/htseq-count --format=bam --mode=intersection-strict --stranded=no %s/%s/%s_sort_byname.bam %s > %s/%s/%s.counts \n" % (args.dump_dir, ext_item_name, item, gtf, args.dump_dir, ext_item_name, item)
	cmd += "rm -rf %s/%s/%s_sort_byname.bam \n" % (args.dump_dir, ext_item_name, item)

# Unmap Count 

if args.unmap_count:
	cmd += "num_unmapped=$(samtools view -c -fox4 %s/%s/%s) \n" % (args.dump_dir, ext_item_name, item)
	cmd += "num_total=$(samtools view -c %s/%s/%s) \n" % (args.dump_dir, ext_item_name, item)
	cmd += "echo \"${num_unmapped} \t ${num_total}\" >> %s/%s/%s_unmapped_ratio.txt \n"  % (args.dump_dir, ext_item_name, item)

if args.repeat:
	cmd += "mkdir %s/%s/repeat_profile \n" % (args.dump_dir, ext_item_name)
	cmd += "python /u/home/h/harryyan/project-eeskin/repeat/rprofile/rprofile.py %s/%s/%s %s/%s/repeat_profile/%s_repeat. \n" % (args.dump_dir, ext_item_name, itemargs.dump_dir, ext_item_name, ext_item_name)


if args.gene: 
	with open(args.gene_list_file) as genes:
		for gene in genes: 
			chromosome = str(gene.split(',')[0])
			name = str(gene.split(',')[3])
			pos_one = str(gene.split(',')[4])
			pos_two = str(gene.split(',')[5])
			cmd += "samtools view -bh %s/%s/%s %s:%s-%s > %s/%s/%s_small_%s.bam \n" % (args.dump_dir, ext_item_name, item, chromosome, pos_one, pos_two, args.dump_dir, ext_item_name, ext_item_name, item)
	cmd += "tar -cvf %s_small_bams.tar *small*.bam --remove-files \n " % (ext_item_name)

if args.immune: 
	with open(args.immune_gene_list_file) as genes:
		for gene in genes: 
			chromosome = str(gene.split(',')[0])
			name = str(gene.split(',')[3])
			pos_one = str(gene.split(',')[4])
			pos_two = str(gene.split(',')[5])
			cmd += "samtools view -bh %s/%s/%s %s:%s-%s > %s/%s/%s_immune_%s.bam \n" % (args.dump_dir, ext_item_name, item, chromosome, pos_one, pos_two, args.dump_dir, ext_item_name, ext_item_name, item)

if args.MT:
	cmd += "samtools view -b %s/%s/%s MT > %s/%s/%s_MT.bam \n" % (args.dump_dir, ext_item_name, item, args.dump_dir, ext_item_name, ext_item_name)

