if [ $# -lt 5 ]
then
echo "GTEx Data Download Script, writteh by Harry Yang - harry2416@gmail.com"
echo "[1] the sample file - e.g. gtex_9777_samples.sample"
echo "[2] the target directory for script"
echo "[3] the download token from NCBI dbGaP - e.g. AF54A1549..."
echo "[4] aspera dir - e.g. /u/home/h/harryyan/.aspera/connect/"
echo "[5] dump dir - the downloaded files will be stored here"
exit 1
fi

pwd
gitdir=/u/home/s/serghei/collab/gtex/bam_processing/
target=$5
mkdir $target
dir=$2
mkdir $dir
cd $dir

while read line
do
item=$(echo $line | awk -F '/' '{print $6}' | awk -F '.ncbi_enc' '{print $1}')
itemname=$(echo $item | awk -F '.' '{print $1}')
echo $line
echo $item
echo $itemname

echo "\"${4}/bin/ascp\" -QTr -l 300M -k 1 -i \"${4}/etc/asperaweb_id_dsa.openssh\" -W ${3} dbtest@gap-upload.ncbi.nlm.nih.gov:data/instant/eeskin/49848/${line} ${target}">> download_${itemname}.sh
echo "$gitdir/bam_process_indiv.sh $5 ${item}.ncbi_enc /u/home/s/serghei/collab/gtex/decryption/ /u/home/s/serghei/collab/gtex/processed/ /u/home/s/serghei/collab/gtex/gtex_ig_C_TCR_C_genes_list.txt  /u/home/s/serghei/gtex/gtex_immune_gene_list.txt" >> download_${itemname}.sh
echo "echo \"DONE ${itemname}\"">> download_${itemname}.sh
done<$1
