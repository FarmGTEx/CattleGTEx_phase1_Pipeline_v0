#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal              # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 24                       # Number of CPU cores
#SBATCH --mem=51200              # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J rna_analysis              # Name of the job
#SBATCH --output=slurm_%x_%A.out   # STDOUT
#SBATCH --error=slurm_%x_%A.out    # STDERR
#SBATCH -t 12:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
# Create a temporary directory for the job in local storage - DO NOT CHANGE #
TMPDIR=/scratch/$USER/$SLURM_JOBID
export TMPDIR
mkdir -p $TMPDIR
#=========================================================================#
#         Your job script                                                 #
#=========================================================================#
# Replace the following with your work to be executed.
echo "Job start at $(date '+%d_%m_%y_%H_%M_%S')"
source "/home/huicongz/miniconda3/etc/profile.d/conda.sh"
      conda activate nf-farmgtex
      
####extract sample and genotype file
sample="REPLACE" 
#sample="Immune System"        
sample=$(echo "$sample" | sed 's/ /_/g')
mkdir -p $sample/discovery/
mkdir -p $sample/validation/ 
Rscript /faststorage/project/farmgtex/pipeline/internal_validation/prepare_annot.r $sample
cd /faststorage/project/farmgtex/pipeline/internal_validation/$sample/discovery/ 
Rscript /faststorage/project/farmgtex/pipeline/internal_validation/prepare.r $sample discovery
cd /faststorage/project/farmgtex/pipeline/internal_validation/$sample/validation/
Rscript /faststorage/project/farmgtex/pipeline/internal_validation/prepare.r $sample validation


###cis eqtl mapping
mkdir -p /faststorage/project/farmgtex/pipeline/internal_validation/result/${sample}/validation/
omiga --mode cis --genotype /faststorage/project/farmgtex/QTL_result/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/internal_validation/${sample}/validation/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/pipeline/internal_validation/result/${sample}/validation/ --verbose --debug --dprop-pc-covar 0.001 --permutations 1000 --multiple-testing clipper
mkdir -p /faststorage/project/farmgtex/pipeline/internal_validation/result/${sample}/discovery/
omiga --mode cis --genotype /faststorage/project/farmgtex/QTL_result/genotype/${sample}/${sample} --phenotype /faststorage/project/farmgtex/pipeline/internal_validation/${sample}/discovery/bed/${sample}.expr_tmm_inv.bed.gz --prefix ${sample} --output-dir /faststorage/project/farmgtex/pipeline/internal_validation/result/${sample}/discovery/ --verbose --debug --dprop-pc-covar 0.001 --permutations 1000 --multiple-testing clipper


echo "Job completed at $(date '+%d_%m_%y_%H_%M_%S')"


#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
