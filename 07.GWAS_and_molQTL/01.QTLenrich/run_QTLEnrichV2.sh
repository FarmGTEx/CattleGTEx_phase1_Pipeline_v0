#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 24                       # Number of CPU cores
#SBATCH --mem=102400                # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -J ieQTL                # Name of the job
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

cd /faststorage/project/cattle_gtexs/script/QTLenrich/QTLEnrich-master/src
main_dir="/faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/state_ieQTL/result/"
for trait_dir in "$main_dir"/*/; do
  trait_name="$(basename "$trait_dir")"
  python3 QTLEnrichV2.py --gwas_file /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/GWAS/${trait_name}.txt \
                             --trait_name ${trait_name} \
                             --qtl_directory /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/tissue_shared/Group6/ \
                             --file_name .cis_qtl_pairs.significant.txt \
                             --qtl_type best_eqtl \
                             --confounders_table /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/tissue_shared/Group1_output_confounder.txt \
                             --null_table /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/tissue_shared/Group1_null_variants_table.txt \
                             --gencode_file /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/Bos_taurus.ARS-UCD1.2.gtf \
                             --gwas_p_value 0.05 \
                             --exp_label eQTL_${trait_name} \
                             --output_directory /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/tissue_shared/result/Group6/${trait_name}/ \
                             --GeneEnrich_input
done
#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID