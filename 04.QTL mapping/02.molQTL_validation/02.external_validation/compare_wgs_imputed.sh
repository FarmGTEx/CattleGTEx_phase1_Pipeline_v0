#!/bin/bash
set -e # 该命令确保如果任何步骤失败，脚本将立即停止

# --- 1. 用户配置区 (您的路径) ---
#!/bin/bash

#!/bin/bash

# --- 1. 定义输入文件和参数 ---
VCF_RNA="/faststorage/project/farmgtex/pipeline/WGS_validate/RNA_temp/recal/recal.final_batch1.bcf"
VCF_ALL_LOCI="/faststorage/project/farmgtex/breeding_project/WGS_vcf/final_vcf/cohort.merged.vcf.gz"
REF_FASTA="/faststorage/project/farmgtex/gtex/cattle_genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
BIM_FILE="/faststorage/project/farmgtex/gtex/VCFfile/chr_final/recal/final_vcf/VCF_file/INFO_0.75_all.bim"
THREADS=8

# --- 2. 定义输出文件名 ---
PROCESSED_RNA="RNA.processed.vcf.gz"
PROCESSED_ALL_LOCI="all_chr_loci.processed.vcf.gz"
WGS_SAMPLES="wgs_samples.txt"
RNA_SAMPLES="rna_samples.txt"
COMMON_SAMPLES="common_samples.txt"
RENAME_MAP="wgs_rename_map.txt"
RENAMED_ALL_LOCI="all_loci.renamed.vcf.gz"
MERGED_VCF="merged_2.vcf.gz"
MERGED_NOMISS_VCF="merged_nomiss.vcf.gz"
FINAL_MERGED_VCF="merged.final.vcf.gz"
SAMPLE_PAIRS_FILE="sample_pairs.txt"
SNP_POS_FILE="snps_to_keep.positions_2.txt"

# --- 3. 脚本主流程 ---
echo "======== [阶段一] 开始预处理和合并VCF文件 (已优化效率) ========"

echo "--> 正在从BIM文件准备SNP位置列表..."
#awk '{print $1"\t"$4}' "$BIM_FILE" > "$SNP_POS_FILE"
bcftools query -f '%CHROM\t%POS\n' $VCF_RNA > "$SNP_POS_FILE"
bcftools query -f '%CHROM\t%POS\t%INFO/INFO\n' $VCF_RNA > new_info.txt
awk '{print $1"_"$2 "\t" $3}' new_info.txt > new_info_merged.txt
echo "    - SNP位置文件准备完成 -> $SNP_POS_FILE"

echo "--> 正在处理 RNA VCF (优化流程)..."
bcftools view -h "$VCF_RNA" | sed 's/_bqsr//g' > rna_new_header.txt
# 优化：先用 view 过滤出需要的SNP位点，再进行 norm
bcftools view -T "$SNP_POS_FILE" --types snps "$VCF_RNA" | \
bcftools norm -f "$REF_FASTA" -c s -m- --threads "$THREADS" - | \
sed 's/|/\//g; s#1/0#0/1#g' | \
bcftools sort  - | \
grep -v '^#' | \
cat rna_new_header.txt - | \
bcftools norm -d exact --threads "$THREADS" - | \
bgzip -@ "$THREADS" > "$PROCESSED_RNA"

# 3. 为已排序和压缩的文件建立索引
echo "    - 正在为输出文件建立索引..."
tabix -p vcf "$PROCESSED_RNA"

echo "    - 处理完成 -> $PROCESSED_RNA"

echo "--> 正在处理 WGS VCF (优化流程)..."
# 优化：先用 view 过滤出需要的SNP位点，再进行后续处理
bcftools view -T "$SNP_POS_FILE" --types snps "$VCF_ALL_LOCI" | \
bcftools annotate -x ^FORMAT/GT | \
bcftools norm -f "$REF_FASTA" -c s -m- --threads "$THREADS" - | \
sed 's/|/\//g; s#1/0#0/1#g' | \
bcftools sort  - | \
bcftools norm -d exact --threads "$THREADS" - | \
bgzip -@ "$THREADS" > "$PROCESSED_ALL_LOCI"


tabix -p vcf "$PROCESSED_ALL_LOCI"
echo "    - 处理完成 -> $PROCESSED_ALL_LOCI"

# (后续所有步骤保持不变)
echo "--> 正在重命名样本并合并文件..."
bcftools query -l "$PROCESSED_ALL_LOCI" > "$WGS_SAMPLES"
awk '{print $1"\t"$1"_wgs"}' "$WGS_SAMPLES" > "$RENAME_MAP"
bcftools reheader -s "$RENAME_MAP" -o "$RENAMED_ALL_LOCI" "$PROCESSED_ALL_LOCI"
tabix -p vcf "$RENAMED_ALL_LOCI"
bcftools merge --threads "$THREADS" -O z -o "$MERGED_VCF" "$PROCESSED_RNA" "$RENAMED_ALL_LOCI"
tabix -p vcf "$MERGED_VCF"
echo "    - 合并完成 -> $MERGED_VCF"

bcftools view -g ^miss "$MERGED_VCF" -Oz -o "$MERGED_NOMISS_VCF"

echo "--> 对合并文件进行最终标准化以拆分新产生的多等位基因位点..."
bcftools norm -m- --threads "$THREADS" -O z -o "$FINAL_MERGED_VCF" "$MERGED_NOMISS_VCF"
tabix -p vcf "$FINAL_MERGED_VCF"
bcftools annotate -x INFO "$FINAL_MERGED_VCF" -Oz -o "${FINAL_MERGED_VCF%.vcf.gz}.noINFO.vcf.gz"
tabix -p vcf "${FINAL_MERGED_VCF%.vcf.gz}.noINFO.vcf.gz"

echo "    - 最终标准化完成 -> $FINAL_MERGED_VCF"

echo "======== [阶段一] 全部处理完成! ========"


# --- 步骤 E: 创建样本配对文件 ---
echo "--> 正在为共同样本创建配对文件..."
bcftools query -l "$PROCESSED_RNA" > "$RNA_SAMPLES"
comm -12 <(sort "$RNA_SAMPLES") <(sort "$WGS_SAMPLES") > "$COMMON_SAMPLES"
awk '{print $1"\t"$1"_wgs"}' "$COMMON_SAMPLES" > "$SAMPLE_PAIRS_FILE"
echo "    - 找到 $(wc -l < $COMMON_SAMPLES) 个共同样本，配对文件创建完成: $SAMPLE_PAIRS_FILE"
echo ""
echo "=====> Linux端流程执行完毕！ <====="
echo "最终产出文件 (供R使用):"
echo "1. VCF文件: $FINAL_MERGED_VCF"
echo "2. 样本配对文件: $SAMPLE_PAIRS_FILE"

# 清理大部分中间文件，但保留处理过的单个VCF用于检查
rm -f rna_new_header.txt wgs_samples.txt rna_samples.txt common_samples.txt wgs_rename_map.txt
rm -f "$RENAMED_ALL_LOCI"* "$MERGED_VCF"*


bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' merged.final.noINFO.vcf.gz > variants_info.txt
bcftools query -f '[%GT\t]\n' merged.final.noINFO.vcf.gz > gt_matrix.txt
bcftools query -l merged.final.noINFO.vcf.gz > samples.txt

