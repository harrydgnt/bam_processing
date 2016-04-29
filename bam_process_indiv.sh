if [ $# -lt 5 ] 
then 
echo "GTEx data process - Written by Harry Yang - harry2416@gmail.com"
echo "[1] dir with .ncbi_enc files" 
echo "[2] sample files for ncbi_enc files"
echo "[3] decryption depository - provided as you set up the decryption key"
echo "	e.g. /u/home/h/harryyan/project-eeskin/decryption_test/"
echo "[4] the dump directory where decrypted bam files would be stored." 
echo "	the dump dir will be created and scripts will be written in this folder!" 
echo "[5] gene list - i.e. gene_coordinate.list"

exit 1
fi 
# $1 = $PWD - download dir 
# $2 = *.ncbi_enc 
# $3 = decryption dir 
# $4 = dump dir - save to this dir 
# $5 = gene list 

pwd 
decdir=$3
dumpdir=$4
mkdir $dumpdir
cd $dumpdir 
pwd 

line=$2
#echo $line 
item=$(echo $line | awk -F '.ncbi_enc' '{print $1}')
itemname=$(echo $item | awk -F '.' '{print $1}') 
ext_item_name=$(echo $item | awk -F '.bam' '{print $1}')
echo $item 

### each sample has dedicated folder

mkdir $ext_item_name
cd $ext_item_name

### decrypt .ncbi_enc file to bam
echo "mv ${1}/${line} ${decdir}" >> run_${itemname}.sh
echo "cd ${decdir}" >> run_${itemname}.sh
echo "${decdir}/vdb-decrypt -f ${decdir}/${line}" >> run_${itemname}.sh 
echo "mv ${decdir}/$item ${dumpdir}/$ext_item_name" >> run_${itemname}.sh 



### HTSeq - default gtf is set

gtf=/u/home/s/serghei/project/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf 

echo "#!/bin/bash" >> run_${itemname}.sh
echo ". /u/local/Modules/default/init/modules.sh" >> run_${itemname}.sh
echo "module load samtools" >> run_${itemname}.sh
echo "module load bamtools" >> run_${itemname}.sh
echo "cd ${dumpdir}/${ext_item_name}" >> run_${itemname}.sh 
echo "samtools sort -n ${dumpdir}/${ext_item_name}/$item ${dumpdir}/${ext_item_name}/${item}_sort_byname" >> run_${itemname}.sh
echo "module load python" >> run_${itemname}.sh
echo "python /u/home/h/harryyan/project-eeskin/utilities/HTSeq-0.6.1/scripts/htseq-count --format=bam --mode=intersection-strict --stranded=no ${dumpdir}/${ext_item_name}/${item}_sort_byname.bam $gtf >${dumpdir}/${ext_item_name}/${item}.counts" >> run_${itemname}.sh
echo "rm -rf ${dumpdir}/${ext_item_name}/${item}_sort_byname.bam" >> run_${itemname}.sh




### bam index

echo "INDEXING THE BAM FILE: ${item}" >> run_${itemname}.sh

echo "samtools index ${dumpdir}/${ext_item_name}/${item}" >> run_${itemname}.sh




### unmapped read

echo "num_unmapped=\$(samtools view -c -fox4 ${dumpdir}/${ext_item_name}/${item})" >> run_${itemname}.sh
echo "num_total=\$(samtools view -c ${dumpdir}/${ext_item_name}/${item})" >> run_${itemname}.sh
echo "echo "\${num_unmapped}	\${num_total}" >> ${dumpdir}/unmapped_ratio.txt" >> run_${itemname}.sh


### repeat profile

echo "mkdir ${dumpdir}/${ext_item_name}/repeat_profile" >> run_${itemname}.sh
echo "python /u/home/h/harryyan/project-eeskin/repeat/rprofile/rprofile.py ${dumpdir}/${ext_item_name}/${item} ${dumpdir}/${ext_item_name}/repeat_profile/${ext_item_name}_repeat." >> run_${itemname}.sh





### bam file extractor for each gene
while read gene
do
chr=$(echo $gene | awk -F ',' '{print $1}')
name=$(echo $gene | awk -F ',' '{print $4}')
pos_one=$(echo $gene | awk -F ',' '{print $5}')
pos_two=$(echo $gene | awk -F ',' '{print $6}')

echo "samtools view -bh ${dumpdir}/${ext_item_name}/${item} $chr:$pos_one-$pos_two > ${dumpdir}/${ext_item_name}/${ext_item_name}_small_${name}.bam" >> run_${itemname}.sh
echo "samtools view -c ${dumpdir}/${ext_item_name}/${ext_item_name}_$name.bam | echo "${name} \$1">>${dumpdir}/${ext_item_name}/${ext_item_name}.count_per_gene" >> run_${itemname}.sh
done<$5

### TODO - ZIP the bam files 
echo "cd ${dumpdir}/${ext_item_name}/" >> run_${itemname}.sh
echo "tar -cvf ${ext_item_name}_small_bams.tar *small*.bam --remove-files" >> run_${itemname}.sh

### extract unmapped bam - unmappped fastq extraction has been disabled
echo "samtools view -b -f 4 ${dumpdir}/${ext_item_name}/${item} > ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.bam" >> run_${itemname}.sh
echo "samtools index ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.bam" >> run_${itemname}.sh
#echo "bamtools convert -in ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.bam -format fastq > ${dumpdir}/${ext_item_name}/${ext_item_name}.fastq">>run_${itemname}.sh
#echo "gzip ${dumpdir}/${ext_item_name}/${ext_item_name}.fastq">>run_${itemname}.sh
#echo "rm ${dumpdir}/${ext_item_name}/${ext_item_name}.bam">>run_${itemname}.sh
#echo "rm ${dumpdir}/${ext_item_name}/${ext_item_name}.bam.bai">>run_${itemname}.sh



### QC STEPS 
#echo "module load bedtools" >> run_${itemname}.sh
#echo "bamToFastq -i ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.bam -fq ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.fastq" >> run_${itemname}.sh 
#echo "/u/home/s/serghei/code/bin/fastq_quality_filter -v -Q 33 -q 20 -p 75 -i ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.fastq -o ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.fastq.trim">>run_${itemname}.sh
#echo "mv ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.fastq.trim ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped_trim.fastq" >> run_${itemname}.sh 
#echo "awk '{if(NR%4==1) {printf(\">%s\n\",substr(\$0,2));} else if(NR%4==2) print;}' ${dumpdir}/${ext_item_name}/${ext_item_name}_trim.fastq >${dumpdir}/${ext_item_name}/${ext_item_name}_trim.fa" >> run_${itemname}.sh 
#echo "/u/home/s/serghei/project/code/import/seqclean-x86_64/seqclean ${dumpdir}/${ext_item_name}/${ext_item_name}_trim.fa -l 50 -M " >> run_${itemname}.sh
#echo "mv ${dumpdir}/${ext_item_name}/${ext_item_name}_trim.fa.clean ${dumpdir}/${ext_item_name}/${ext_item_name}_trim_clean.fa" >> run_${itemname}.sh 
#echo "echo \"$ext_item_name DONE \" > ${dumpdir}/${ext_item_name}/${ext_item_name}_QC.DONE" >> run_${itemname}.sh  
#echo "rm ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.bam" >> run_${itemname}.sh 
#echo "rm ${dumpdir}/${ext_item_name}/${ext_item_name}_trim.fastq" >> run_${itemname}.sh 
#echo "echo \"${ext_item_name}.unmap_process.DONE\">> ${dumpdir}/${ext_item_name}/${ext_item_name}_unmap.DONE" >>run_${itemname}.sh
#echo "rm ${dumpdir}/${ext_item_name}/$item">>run_${itemname}.sh

### add ROP 

echo "python /u/home/h/harryyan/project-eeskin/rop/rop.py --b --qsubArray ${dumpdir}/${ext_item_name}/${ext_item_name}.unmapped.bam ${dumpdir}/${ext_item_name}/rop/" >> run_${itemname}.sh

### TODO - remove intermediate files


### genomic profile
echo "mkdir ${dumpdir}/${ext_item_name}/genomic_profile">>run_${itemname}.sh
echo "python /u/home/h/harryyan/project-eeskin/repeat/gprofile/gprofilePE.py --readPerCategory ${dumpdir}/${ext_item_name}/${item} ${dumpdir}/${ext_item_name}/genomic_profile/${ext_item_name}_genome h">>run_${itemname}.sh
echo "rm ${dumpdir}/${ext_item_name}/${ext_item_name}.bam">>run_${itemname}.sh
echo "rm ${dumpdir}/${ext_item_name}/${ext_item_name}.bam.bai">>run_${itemname}.sh
echo "echo \"${ext_item_name}.gprofile.DONE\">> ${dumpdir}/${ext_item_name}/${ext_item_name}_gprofile.DONE" >>run_${itemname}.sh

qsub -cwd -V -N run_${itemname} -l h_data=8G,time=24:00:00 run_${itemname}.sh 

cd ..
