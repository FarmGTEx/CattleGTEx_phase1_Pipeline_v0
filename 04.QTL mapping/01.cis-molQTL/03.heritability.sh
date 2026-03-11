#!/bin/bash

echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
      
####extract sample and genotype file
sample=$1     
sample=$(echo "$sample" | sed 's/ /_/g')

### cis-eQTL heri
mkdir -p /faststorage/project/farmgtex/QTL_result_new/cis_her/${sample}/
omiga --mode her_est --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/eQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample}_2 --output-dir /faststorage/project/farmgtex/QTL_result_new/cis_her/${sample}/ --h2-model Ac --prefix ${sample} --dprop-pc-covar 0.001 --threads 96

mkdir -p /faststorage/project/farmgtex/QTL_result_new/her/cis/eeQTL/
omiga --mode her_est --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/eeQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/her/cis/eeQTL/${sample}/ --h2-model Ac --prefix ${sample} --dprop-pc-covar 0.001 --threads 24 --pheno-group /faststorage/project/farmgtex/pipeline/eeQTL/${sample}/bed/${sample}.group

mkdir -p /faststorage/project/farmgtex/QTL_result_new/her/cis/eeQTL/
omiga --mode her_est --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/enQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/her/cis/enQTL/${sample}/ --h2-model Ac --prefix ${sample} --dprop-pc-covar 0.001 --threads 12

mkdir -p /faststorage/project/farmgtex/QTL_result_new/her/cis/eeQTL/
omiga --mode her_est --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/isoQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/her/cis/isoQTL/${sample}/ --h2-model Ac --prefix ${sample} --dprop-pc-covar 0.001 --threads 12

mkdir -p /faststorage/project/farmgtex/QTL_result_new/her/cis/eeQTL/
omiga --mode her_est --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype --phenotype /faststorage/project/farmgtex/pipeline/sQTL/${sample}/bed/${sample}_filted_qqnorm.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/her/cis/sQTL/${sample}/ --h2-model Ac --prefix ${sample} --dprop-pc-covar 0.001 --threads 16 

mkdir -p /faststorage/project/farmgtex/QTL_result_new/her/cis/eeQTL/
omiga --mode her_est --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/stQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/QTL_result_new/her/cis/stQTL/${sample}/ --h2-model Ac --prefix ${sample} --dprop-pc-covar 0.001 --threads 12

mkdir -p /faststorage/project/farmgtex/QTL_result_new/her/cis/eeQTL/
omiga --mode her_est --genotype /faststorage/project/farmgtex/QTL_result_new/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/3aQTL/${sample}/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir  /faststorage/project/farmgtex/QTL_result_new/her/cis/3aQTL/${sample} --h2-model Ac --prefix ${sample} --dprop-pc-covar 0.001 --threads 12