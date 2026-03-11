#!/bin/bash

number=10000
selector="PLSR"
method="svp_p"

cd ${work_path} && mkdir -p predict_result && cd predict_result

#Pick informative SNPs for breed prediction
Rscript ${work_path}/pick_snp.R $number

#Pick samples for breed prediction
Rscript ${work_path}/Predict_Sample_Extract.R 

plink --bfile ${work_path}/predict_samples  --recode A --extract ${selector}_${number}snp.txt --threads 12 --cow --out ${selector}_${number}_ped
python ${work_path}/Predict.py ${number} ${selector} ${method} 