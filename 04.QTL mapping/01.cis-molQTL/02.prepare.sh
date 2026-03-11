#!/bin/bash

##eQTL
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
sample=$1
sample=$(echo "$sample" | sed 's/ /_/g')
mkdir -p $sample
Rscript prepare_annot.r $sample
cd $sample
cat sample_list.txt | awk '{if ($0 ~ /-/) print $0 ".bqsr"; else print $0 "_bqsr"}' > VCF_keep.txt
awk '{print $1, $1}' VCF_keep.txt > VCF_keep_plink.txt
Rscript /faststorage/project/farmgtex/pipeline/eQTL/prepare_eQTL.r $sample
Rscript /faststorage/project/farmgtex/pipeline/eeQTL/prepare_eQTL_annot.r $sample

mkdir -p /faststorage/project/farmgtex/QTL_result_new/genotype/$sample
cd /faststorage/project/farmgtex/QTL_result_new/genotype/$sample
bcftools view -S /faststorage/project/farmgtex/pipeline/eQTL/${sample}/VCF_keep.txt  /faststorage/project/farmgtex/gtex/VCFfile/chr_final/recal/final_vcf/VCF_file/INFO_0.75_all.vcf -o /faststorage/project/farmgtex/QTL_result/genotype/$sample/genotype.vcf
plink --vcf /faststorage/project/farmgtex/QTL_result/genotype/$sample/genotype.vcf --make-bed --maf 0.05 --mac 6 --out /faststorage/project/farmgtex/QTL_result/genotype/$sample/$sample --cow --double-id
awk 'BEGIN {OFS="\t"} { $2 = $1"_"$4; print $0 }' ${sample}.bim > temp.bim && mv temp.bim ${sample}.bim
plink \
  --bfile /faststorage/project/farmgtex/gtex/VCFfile/chr_final/recal/final_vcf/VCF_file/INFO_0.75_all \
  --keep /faststorage/project/farmgtex/pipeline/eQTL/${sample}/VCF_keep_plink.txt\
  --make-bed --keep-allele-order \
  --cow \
  --out /faststorage/project/farmgtex/QTL_result_new/genotype/$sample/$sample
sed -i 's/_bqsr//g; s/.bqsr//g' ${sample}.fam

##eeQTL
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
sample=$1       
sample=$(echo "$sample" | sed 's/ /_/g')
mkdir -p $sample
cd $sample
Rscript /faststorage/project/farmgtex/pipeline/eeQTL/prepare_eeQTL_annot.r $sample
Rscript /faststorage/project/farmgtex/pipeline/eeQTL/prepare_eeQTL.r $sample
rm /faststorage/project/farmgtex/pipeline/eeQTL/${sample}/*rds

##enQTL
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
sample=$1        
sample=$(echo "$sample" | sed 's/ /_/g')
mkdir -p $sample
cd $sample
Rscript /faststorage/project/farmgtex/pipeline/enQTL/prepare_enQTL_annot.r $sample
Rscript /faststorage/project/farmgtex/pipeline/enQTL/prepare_enQTL.r $sample
rm /faststorage/project/farmgtex/pipeline/enQTL/${sample}/*rds

##isoQTL
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
sample=$1      
sample=$(echo "$sample" | sed 's/ /_/g')
mkdir -p $sample
Rscript prepare_annot.r $sample
cd $sample
Rscript /faststorage/project/farmgtex/pipeline/isoQTL/prepare_isoQTL.r $sample
Rscript /faststorage/project/farmgtex/pipeline/isoQTL/prepare_isoQTL_annot.r $sample

##sQTL
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
sample=$1        
sample=$(echo "$sample" | sed 's/ /_/g')
mkdir -p $sample
cd $sample
Rscript /faststorage/project/farmgtex/pipeline/sQTL/prepare_sQTL_annot.r $sample
Rscript /faststorage/project/farmgtex/pipeline/sQTL/prepare_sQTL.r $sample
rm /faststorage/project/farmgtex/pipeline/sQTL/${sample}/*rds

##stQTL
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
sample=$1        
sample=$(echo "$sample" | sed 's/ /_/g')
mkdir -p $sample
cd $sample
Rscript /faststorage/project/farmgtex/pipeline/stQTL/prepare_stQTL_annot.r $sample
Rscript /faststorage/project/farmgtex/pipeline/stQTL/prepare_stQTL.r $sample
rm /faststorage/project/farmgtex/pipeline/stQTL/${sample}/*rds

##3a'QTL
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
sample=$1        
sample=$(echo "$sample" | sed 's/ /_/g')
mkdir -p $sample
cd $sample
Rscript /faststorage/project/farmgtex/pipeline/3aQTL/prepare_3aQTL_annot.r $sample
Rscript /faststorage/project/farmgtex/pipeline/3aQTL/prepare_3aQTL.r $sample
rm /faststorage/project/farmgtex/pipeline/3aQTL/${sample}/*rds
