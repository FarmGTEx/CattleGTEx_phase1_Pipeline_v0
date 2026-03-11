#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#

main_dir="/faststorage/project/farmgtex/Population_eQTL/eQTL_dairy_beef/genotype/"
geno_dir="/faststorage/project/farmgtex/Population_eQTL/eQTL_dairy_beef/genotype/"

for subdir in "$main_dir"/*/; do
  cd "${geno_dir}/$(basename "$subdir")"
  bfile_prefix="${geno_dir}/$(basename "$subdir")/$(basename "$subdir")_bfile_admix.bed"
  for k in {2..3} ; do
    sbatch /faststorage/project/cattle_gtexs/script/admixture.sh "$bfile_prefix" "$k"
  done
done

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID