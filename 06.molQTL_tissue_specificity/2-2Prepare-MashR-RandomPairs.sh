#!/bin/bash
dir_nominal=${1}
dir_output=${2}
tis_i=${3}
CHR=${4}

mkdir -p ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 2. extract all nominal pairs from nominal results for each tissue
# permutation results
filter_list=$(awk -F'/' '{print $NF}' "$dir_output/permutation_files_filter.txt")
perm_files=$(awk -F'/' '{print $NF}' "$dir_output/permutation_files_filter.txt")
NAMEs=()
while IFS= read -r line; do
    name=$(basename "$line" ".cis_qtl.txt.gz")
    NAMEs+=("$name")
done < <(grep "${dir_output}/permutation_files_filter.txt" -f "/faststorage/project/farmgtex/pipeline/eQTL_tissues" )

name=${NAMEs[tis_i]}
nominal_files=($(find "${dir_nominal}" \
    -name "${name}*.cis_qtl_pairs.${CHR}.txt.gz" \
    ! -name "*LMM*"))
rm -f ${dir_output}/${name}_Chr${CHR}.nominal_files.txt
for l in ${nominal_files[*]}
do
{
    echo ${l} >> ${dir_output}/${name}_Chr${CHR}.nominal_files.txt
}
done
# extract_pairs
python3 /faststorage/project/farmgtex/pipeline/code/mashR//extract_pairs_tjy.py ${dir_output}/${name}_Chr${CHR}.nominal_files.txt ${dir_output}/nominal_pairs.combined_signifpairs.txt.gz ${name}_Chr${CHR}.nominal_pairs --chrom ${CHR} -o ${dir_output}
# > output file: *_nominal_pairs.extracted_pairs.txt.gz
rm -f ${dir_output}/${name}_Chr${CHR}.nominal_files.txt
