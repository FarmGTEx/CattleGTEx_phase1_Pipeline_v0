#!/bin/bash
#--------------------------------------------------------------------------#
#              Edit Job specifications                                     #
#--------------------------------------------------------------------------#
#SBATCH -p normal                 # Name of the queue
#SBATCH -N 1                       # Number of nodes(DO NOT CHANGE)
#SBATCH -n 12                       # Number of CPU cores
#SBATCH --mem=512000               # Memory in MiB(10 GiB = 10 * 1024 MiB)
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

cd /faststorage/project/cattle_gtexs/script/QTLenrich/

python3 extract_unique_variants.py --directory /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2/ --file_name .cis_qtl_pairs.significant.txt --output_file /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_significant_variants.txt

python3 create_null_variants.py --all_variants_file /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/all_variant.txt --significant_variants_file /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_significant_variants.txt --output_file /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_null_variants.txt

python3 Generate_Null_Table.py --QTL_Directory /faststorage/project/farmgtex/Figure6/QTLenrich/eQTL/all_pairs/ --File_Extension .eqtl_allpairs.txt --Null_Variants_List /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_null_variants.txt --Output_File /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_null_variants_table.txt.gz

python3 Generate_Confounders_Table.py --QTL_Directory /faststorage/project/farmgtex/Figure6/QTLenrich/eQTL/all_pairs/ --File_Extension .eqtl_allpairs.txt --Variants_List /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/all_variant.txt --LD_Proxy /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/number_ld_per_variant.txt --GENCODE_File /faststorage/project/cattle_gtexs/CattleGTEx/OmiGA/QTLenrich/ct_seQTL/Bos_taurus.ARS-UCD1.2.gtf --Output_File /faststorage/project/farmgtex/Figure6/QTLenrich/3aQTL/Group2_output_confounder.txt.gz

#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID