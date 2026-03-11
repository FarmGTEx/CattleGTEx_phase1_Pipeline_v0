#!/bin/bash

input="TPM_median_tissue_matrix.csv"

metrics=("counts" "tau" "gini" "simpson" "shannon_specificity" \
         "roku_specificity" "tsi" "zscore" "spm" "spm_dpm" \
         "js_specificity" "js_specificity_dpm")

for method in "${metrics[@]}"; do
    echo "Running tspex for metric: $method"
    tspex "$input" "tissue_${method}.txt" "$method"
done

echo "All tspex metrics finished!"
