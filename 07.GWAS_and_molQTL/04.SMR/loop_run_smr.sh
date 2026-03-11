#!/bin/bash

main_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/state_ieQTL/result/"
for trait_dir in "$main_dir"/*/; do
  trait="$(basename "$trait_dir")"
  sbatch /faststorage/project/farmgtex/GWAS_coloc/SMR/run_smr.sh "$trait"
done
