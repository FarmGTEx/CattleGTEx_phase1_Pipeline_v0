

bcftools view -i 'INFO>0.75 && TYPE="snp"' /faststorage/project/farmgtex/pipeline/QTL_validation/recal.final_batch1.bcf \
  -Oz -o RNA.vcf.gz
  
# /faststorage/project/farmgtex/pipeline/QTL_validation/recal.final_batch1.bcf
#"/faststorage/project/farmgtex/breeding_project/WGS_vcf/final_vcf/cohort.merged.vcf.gz"  

bcftools query -l /faststorage/project/farmgtex/breeding_project/WGS_vcf/final_vcf/cohort.merged.vcf.gz \
 > wgs_samples.txt
  
bcftools query -l /faststorage/project/farmgtex/pipeline/QTL_validation/recal.final_batch1.bcf \
  | sed 's/_bqsr$//' \
  > rna_samples.txt
  
comm -12 <(sort wgs_samples.txt) <(sort rna_samples.txt) > overlap_samples.txt  
  
bcftools view -v snps /faststorage/project/farmgtex/breeding_project/WGS_vcf/final_vcf/cohort.merged.vcf.gz \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > wgs_snps.txt

# RNA SNP
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' /faststorage/project/farmgtex/pipeline/QTL_validation/RNA.vcf.gz\
  > rna_snps.txt  
  
  
sort wgs_snps.txt > wgs_snps.sorted.txt
sort rna_snps.txt > rna_snps.sorted.txt
comm -12 wgs_snps.sorted.txt rna_snps.sorted.txt > overlap_snps.txt

grep -F -f overlap_samples.txt <(bcftools query -l RNA.vcf.gz) \
  > overlap_rna_samples.txt
  
  
  bcftools view \
  -S overlap_samples.txt \
  -T overlap_snps.txt \
  /faststorage/project/farmgtex/breeding_project/WGS_vcf/final_vcf/cohort.merged.vcf.gz | \
bcftools annotate \
  -x INFO,FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PL,FORMAT/PGT,FORMAT/PID,FORMAT/PS \
  -O z \
  -o WGS.filtered.vcf.gz

  
  bcftools view \
  -S overlap_rna_samples.txt \
  -T overlap_snps.txt \
  -O z \
  -o RNA.filtered.vcf.gz \
  RNA.vcf.gz
bcftools index WGS.filtered.vcf.gz

# 为 RNA VCF 创建索引
bcftools index RNA.filtered.vcf.gz 

bcftools merge \
  WGS.filtered.vcf.gz \
  RNA.filtered.vcf.gz \
  -O z -o merged_filtered.vcf.gz
  
  bcftools view \
  -S overlap_rna_samples.txt \
  -O z \
  -o RNA.filtered_all.vcf.gz \
  RNA.vcf.gz  
  
  
  bcftools view \
  -S overlap_samples.txt \
  /faststorage/project/farmgtex/breeding_project/WGS_vcf/final_vcf/cohort.merged.vcf.gz | \
bcftools annotate \
  -x INFO,FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PL,FORMAT/PGT,FORMAT/PID,FORMAT/PS \
  -O z \
  -o WGS.filtered_all.vcf.gz
 
 
 
 
 bcftools index WGS.filtered_all.vcf.gz

# 为 RNA VCF 创建索引
bcftools index RNA.filtered_all.vcf.gz 

  bcftools merge \
  WGS.filtered_all.vcf.gz \
  RNA.filtered_all.vcf.gz \
  -O z -o merged_filtered_all.vcf.gz
  
plink --vcf merged_filtered_all.vcf.gz --make-bed --out merged  --cow
  