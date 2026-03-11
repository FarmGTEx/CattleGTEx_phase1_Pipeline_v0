#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal              # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 8                       # Number of CPU cores
#SBATCH --mem=50g              # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J rna_analysis              # Name of the job
#SBATCH --output=slurm_%x_%A.out   # STDOUT
#SBATCH --error=slurm_%x_%A.out    # STDERR
#SBATCH -t 12:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
#SBATCH --account=farmgtex
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

cd /faststorage/project/farmgtex/gtex/bqsr_WGS/

gatk CombineGVCFs -R /faststorage/project/farmgtex/gtex/cattle_genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa $(cat /faststorage/project/farmgtex/pipeline/WGS_validate/gvcf_args.txt) -O /faststorage/project/farmgtex/pipeline/WGS_validate/WGS_data/cohort.g.vcf.gz


echo "Job completed at $(date '+%d_%m_%y_%H_%M_%S')"


#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
