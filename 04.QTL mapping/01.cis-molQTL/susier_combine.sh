#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal              # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 4                       # Number of CPU cores
#SBATCH --mem=20480              # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J susier_conbine_TIS              # Name of the job
#SBATCH --output=PATHsusier_combine_TIS.%J.log   # STDOUT
#SBATCH --error=PATHsusier_combine_TIS.%J.err   # STDERR
#SBATCH -t 4:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR 
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#
# Replace the following with your work to be executed.

source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
conda activate nf-farmgtex
tis=$1
QTL_dir=$2
work_dir=$3


result_dir="${work_dir}/${tis}/"
final_dir="${result_dir}/result/"
QTL_result_dir="${QTL_dir}/"
gene_list_dir="${work_dir}/gene_list/${tis}/"
rm -f ${final_dir}/${tis}.crediblenum.txt ${final_dir}/${tis}.causalnum.txt
ls ${gene_list_dir}/x?? | xargs -n 1 basename | while read num
do
	echo $num
	if [ $(wc -l < ${final_dir}/${num}.susier.gz) -gt 0 ]; then
  # number of credible sets for each gene
  echo ${final_dir}
  ls ${final_dir}/
	zcat ${final_dir}/${num}.susier.gz | awk -v FS="\t" '(NF>=13){print $1"\t"$2"\t"$4}' | sed '1d' | sort -uV | bedtools groupby -i - -g 1,2 -c 3 -o count | awk '{print $1"\t"$2"\t"$3-1}' >> ${final_dir}/${tis}.crediblenum.txt
  # number of causal variants for each credible set
	zcat ${final_dir}/${num}.susier.gz | awk -v FS="\t" '(NF>=13&&$4!=""){print $1"\t"$2"\t"$4}' | sed '1d' | sort | uniq -c >> ${final_dir}/${tis}.causalnum.txt
else
		echo "${tis}/results/susier/${num}.susier.gz is empty"
	fi
done


# combine susier results

zcat ${final_dir}/x??.susier.gz | awk 'NR==1||(NF>=13&&$1!="tissue")' | gzip -c > ${final_dir}/${tis}.susier.gz
# extract 95% credible sets
zcat ${final_dir}/${tis}.susier.gz | awk -v FS="\t" '$4!=""' | gzip -c > ${final_dir}/${tis}.susier.credible.gz
# extract eQTLs in 95% credible sets
csvtk join -t -f 'pheno_id,SNP;pheno_id,variant_id' <(zcat ${final_dir}/${tis}.susier.credible.gz) <(cat ${QTL_result_dir}/${tis}.cis_qtl_pairs.all_significant.txt | cut -f1,2,3,5,7) | gzip -c > ${final_dir}/${tis}.susier.credible.sig.gz
# extract the lead eQTL in each 95% credible set
zcat ${final_dir}/${tis}.susier.credible.sig.gz | awk '{if (NR==1){print}else{key=$2" "$4;if ($5>max[key]){max[key]=$5;line[key]=$0}}}END{for (k in line){print line[k]}}' | gzip -c > ${final_dir}/${tis}.susier.credible.sig.lead.gz

echo "Job completed at $(date '+%d_%m_%y_%H_%M_%S')"


#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
