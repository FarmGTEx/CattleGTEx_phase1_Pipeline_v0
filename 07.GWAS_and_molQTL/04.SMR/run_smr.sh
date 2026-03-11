#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 24                       # Number of CPU cores
#SBATCH --mem=102400                # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J SMR                # Name of the job
#SBATCH --output=slurm_%A.out   # STDOUT
#SBATCH --error=slurm_%A.err    # STDERR
#SBATCH -t 10:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#

set -euo pipefail
trait_name="$1"

main_dir="/faststorage/project/farmgtex/GWAS_coloc/rank_QTL/eQTL/SMR_tissue/"
for tissue_dir in "$main_dir"/*/; do
  tissue_name="$(basename "$tissue_dir")"
  /home/holi/farmgtex/GWAS_coloc/SMR/smr-1.4.0-linux-x86_64/smr --bfile /faststorage/project/farmgtex/QTL_result_new/genotype/${tissue_name}/${tissue_name} --gwas-summary /faststorage/project/farmgtex/GWAS_coloc/SMR/GWAS/${trait_name}_cojo.txt --beqtl-summary /faststorage/project/farmgtex/GWAS_coloc/SMR/BESD_file/${tissue_name}.eqtl_allpairs --out /faststorage/project/farmgtex/GWAS_coloc/SMR/result/${trait_name}/${tissue_name} --thread-num 10 --diff-freq-prop 1
done

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID