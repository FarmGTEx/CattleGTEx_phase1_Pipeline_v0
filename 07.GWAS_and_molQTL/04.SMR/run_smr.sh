#!/bin/bash

set -euo pipefail
trait_name="$1"

main_dir="/faststorage/project/farmgtex/GWAS_coloc/rank_QTL/eQTL/SMR_tissue/"
for tissue_dir in "$main_dir"/*/; do
  tissue_name="$(basename "$tissue_dir")"
  /home/holi/farmgtex/GWAS_coloc/SMR/smr-1.4.0-linux-x86_64/smr --bfile /faststorage/project/farmgtex/QTL_result_new/genotype/${tissue_name}/${tissue_name} --gwas-summary /faststorage/project/farmgtex/GWAS_coloc/SMR/GWAS/${trait_name}_cojo.txt --beqtl-summary /faststorage/project/farmgtex/GWAS_coloc/SMR/BESD_file/${tissue_name}.eqtl_allpairs --out /faststorage/project/farmgtex/GWAS_coloc/SMR/result/${trait_name}/${tissue_name} --thread-num 10 --diff-freq-prop 1
done
