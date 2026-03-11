#!/bin/bash

# Replace the following with your work to be executed.
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
      
####extract sample and genotype file
sample=$1        
sample=$(echo "$sample" | sed 's/ /_/g')

cd /faststorage/project/farmgtex/QTL_result_new/eQTL/${sample}/
omiga --mode cis_independent --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/eQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/eQTL/${sample}/ --verbose --debug --dprop-pc-covar 0.001 --cis-file /faststorage/project/farmgtex/QTL_result_new/eQTL/${sample}/${sample}.cis_qtl.txt.gz --maf-threshold 0.05

cd /faststorage/project/farmgtex/QTL_result_new/eeQTL/${sample}/
omiga --mode cis_independent --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/eeQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/eeQTL/${sample}/ --verbose --debug --dprop-pc-covar 0.001 --cis-file /faststorage/project/farmgtex/QTL_result_new/eeQTL/${sample}/${sample}.cis_qtl.txt.gz

cd /faststorage/project/farmgtex/QTL_result_new/enQTL/${sample}/
omiga --mode cis_independent --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/enQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --cis-window 100000 --output-dir /faststorage/project/farmgtex/QTL_result_new/enQTL/${sample}/ --verbose --debug --dprop-pc-covar 0.001 --cis-file /faststorage/project/farmgtex/QTL_result_new/enQTL/${sample}/${sample}.cis_qtl.txt.gz

cd /faststorage/project/farmgtex/QTL_result_new/isoQTL/${sample}/
omiga --mode cis_independent --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/isoQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/isoQTL/${sample}/ --verbose --debug --dprop-pc-covar 0.001 --cis-file /faststorage/project/farmgtex/QTL_result_new/isoQTL/${sample}/${sample}.cis_qtl.txt.gz

cd /faststorage/project/farmgtex/QTL_result_new/sQTL/${sample}/
omiga --mode cis_independent --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/sQTL/${sample}/bed/${sample}_filted_qqnorm.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/sQTL/${sample}/ --verbose --debug --dprop-pc-covar 0.001 --pheno-group /faststorage/project/farmgtex/pipeline/sQTL/${sample}/sample_group.txt --cis-file /faststorage/project/farmgtex/QTL_result_new/sQTL/${sample}/${sample}.cis_qtl.txt.gz --maf-threshold 0.05

cd /faststorage/project/farmgtex/QTL_result_new/stQTL/${sample}/
omiga --mode cis_independent --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/stQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/stQTL/${sample}/ --verbose --debug --dprop-pc-covar 0.001 --cis-file /faststorage/project/farmgtex/QTL_result_new/stQTL/${sample}/${sample}.cis_qtl.txt.gz 

cd /faststorage/project/farmgtex/QTL_result_new/3aQTL/${sample}/
omiga --mode cis_independent --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/3aQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/3aQTL/${sample}/ --verbose --debug --dprop-pc-covar 0.001 --cis-file /faststorage/project/farmgtex/QTL_result_new/3aQTL/${sample}/${sample}.cis_qtl.txt.gz