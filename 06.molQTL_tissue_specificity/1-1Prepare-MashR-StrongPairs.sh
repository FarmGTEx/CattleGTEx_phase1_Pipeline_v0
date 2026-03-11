#!/bin/bash
dir_nominal=${1}
dir_output=${2}


mkdir -p ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# file list of permutation results
rm -f ${dir_output}/permutation_files.txt
# permutation results
perm_files=($(find "${dir_nominal}" -name "*.cis_qtl.txt.gz" ! -path "*LMM*"))
for l in ${perm_files[*]}
do
{
    echo ${l} >> ${dir_output}/permutation_files.txt
}
done
cat ${dir_output}/permutation_files.txt | grep -f "/faststorage/project/farmgtex/pipeline/eQTL_tissues" > ${dir_output}/permutation_files_filter.txt
python3 /faststorage/project/farmgtex/pipeline/code/mashR/combine_signif_pairs_tjy.py ${dir_output}/permutation_files_filter.txt strong_pairs -o ${dir_output}
#> output file: strong_pairs.combined_signifpairs.txt.gz
