QTL_dir="/faststorage/project/farmgtex/QTL_result/eQTL/"
result_dir="/faststorage/project/farmgtex/pipeline/code/mashR/mashR_new"
QTL_design_file="/faststorage/project/farmgtex/pipeline/QTL_design"

bash 1-1Prepare-MashR-StrongPairs.sh $QTL_dir $result_dir
bash 1-2Prepare-MashR-StrongPairs.sh $QTL_dir $result_dir
bash 1-3Prepare-MashR-StrongPairs.sh $QTL_dir $result_dir

bash 2-1Prepare-MashR-RandomPairs.sh $QTL_dir $result_dir
for tissue in {0..35} 
do 
  for chr in {1..29}  
  do
    #bash 2-2Prepare-MashR-RandomPairs.sh $QTL_dir $result_dir $tissue $chr
    sbatch -p normal -n 1 --mem 30g 2-2Prepare-MashR-RandomPairs.sh $QTL_dir $result_dir $tissue $chr
  done
done

for chr in {1..29}  
  do
    #bash 2-3Prepare-MashR-RandomPairs.sh $result_dir  $chr
    sbatch -p normal -n 1 --mem 30g 2-3Prepare-MashR-RandomPairs.sh  $result_dir $chr
  done
  
conda activate nf-farmgtex
Rscript MashR-random_subset.R $result_dir 1000000

Rscript run_MashR.R  strong_pairs.MashR_input.txt.gz MashR.random_subset.RDS 1 $result_dir 

sbatch -p normal -J run_mashR -n 1 --mem 100g -t 600 --wrap="Rscript run_MashR.R  strong_pairs.MashR_input.txt.gz MashR.random_subset.RDS 0 $result_dir/ " 