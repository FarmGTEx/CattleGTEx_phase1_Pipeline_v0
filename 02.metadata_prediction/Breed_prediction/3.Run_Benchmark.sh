#!/bin/bash

all_count=("50" "100" "200" "500" "1000" "2000" "3000" "5000" "10000" "20000" "30000" "50000")
all_meth=("svp_p" "svp_r" "rf" "xgboost" "knn")
all_sele=("AED" "FST" "PLSR")

for number in "${all_count[@]}"; do
    if ls "${work_path}/${number}"/*.raw 1> /dev/null 2>&1; then
        echo "Raw file exists for $number, skipping..."
    else
	#Extract SNPs and convert the file to PLINK .raw format
        Rscript Pick_SNP.R
    fi	
    for method in "${all_meth[@]}"; do
        for selector in "${all_sele[@]}"; do
            mkdir -p "${work_path}/${number}"
            cd "${work_path}/${number}" || exit 1
            python ML.py "${number}" "${selector}" "${method}"
        done
    done
done 
