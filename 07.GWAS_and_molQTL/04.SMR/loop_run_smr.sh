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

main_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/state_ieQTL/result/"
for trait_dir in "$main_dir"/*/; do
  trait="$(basename "$trait_dir")"
  sbatch /faststorage/project/farmgtex/GWAS_coloc/SMR/run_smr.sh "$trait"
done

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID