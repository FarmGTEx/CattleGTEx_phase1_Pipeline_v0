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
#SBATCH -t 6:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
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
#step3 160g
bin=$1
chr=$2
step=$3
function main() {
        step${step}
}

function step1() {


cd /faststorage/project/farmgtex/gtex/cattle_ref/CattleSNPs
ls *.vcf.gz | while read file
do
cd /faststorage/project/farmgtex/gtex/cattle_ref/CattleSNPs
tabix -p vcf ${file}
bcftools +fill-tags -Oz -o /faststorage/project/farmgtex/gtex/cattle_ref/SNPs/AN_AC.${file} -- $file -t AN,AC
cd /faststorage/project/farmgtex/gtex/cattle_ref/SNPs/
bcftools index AN_AC.${file}
done
}

function step2() {
conda activate glimpse2
cd /faststorage/project/farmgtex/gtex/cattle_ref/SNPs/
cat /faststorage/project/farmgtex/gtex/cattle_ref/chrlist | while read chr
do
GLIMPSE2_chunk --input AN_AC.${chr}.phased.vcf.gz --region ${chr} --window-mb 5.0 --buffer-mb 1.0 --output chunks.${chr}.txt --sequential --threads 12
while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    GLIMPSE2_split_reference --reference AN_AC.${chr}.phased.vcf.gz --input-region ${IRG} --output-region ${ORG} --output ./split/split
done < chunks.$chr.txt
done
cd /faststorage/project/farmgtex/gtex/cattle_ref/SNPs/split/
for bin in `ls -1v split_*.bin` ; do basename $bin .bin ; done | sed 's/^split\_//g' > bins.txt
}


function step3() {
conda activate glimpse2
cd /faststorage/project/farmgtex/pipeline/WGS_validate/RNA_temp2/bqsr_bam/bam_dir/
ls *bam > /faststorage/project/farmgtex/pipeline/WGS_validate/RNA_temp2/bqsr_bam/bam_dir/RNA_list
GLIMPSE2_phase \
            --bam-list /faststorage/project/farmgtex/pipeline/WGS_validate/RNA_temp2/bqsr_bam/bam_dir/RNA_list --reference /faststorage/project/farmgtex/gtex/cattle_ref/SNPs/split/split_${bin}.bin \
            --output /faststorage/project/farmgtex/pipeline/WGS_validate/RNA_temp2/${bin}.bcf --threads 12
}

function step4() {
conda activate glimpse2
cd /faststorage/project/farmgtex/pipeline/WGS_validate/RNA_temp2/
mkdir -p recal
ls -1v  ${chr}_*.bcf >  ${chr}.txt

        GLIMPSE2_ligate \
            --input ${chr}.txt \
            --output impute${chr}.bcf \
            --threads 12
            
conda activate nf-farmgtex            
python3 /faststorage/project/farmgtex/gtex/VCFfile/recalINFO.py --invcf impute${chr}.bcf > ./recal/recal_${chr}.vcf
cd ./recal/
bcftools annotate -x ID,FORMAT/GP,FORMAT/DS,INFO/RAF  -O b -o recal.simplified_${chr}.bcf recal_${chr}.vcf --threads 12
bcftools index recal.simplified_${chr}.bcf --threads 12        
#bcftools view -i 'INFO > 0.6 & AF > 0.05 & AF < 0.95' recal.simplified_${chr}.bcf -Ob -o recal.simplified_filter_${chr}.bcf --threads 12
#bcftools index recal.simplified_filter_${chr}.bcf --threads 12
#bcftools query -f '%CHROM\t%POS\t%INFO\n' recal.simplified_filter_${chr}.bcf > ${chr}.txt
}
function step5() {
conda activate nf-farmgtex
cd /faststorage/project/farmgtex/pipeline/WGS_validate/RNA_temp2/recal/
ls -1v recal.simplified_[0-9]*.bcf >  impute.txt

bcftools index recal.simplified_X.bcf
conda activate glimpse2
          GLIMPSE2_ligate \
          --input impute.txt \
          --output recal.final_batch2.bcf \
          --threads 12
}

function step7() {
conda activate nf-farmgtex 
cd /faststorage/project/farmgtex/pipeline/WGS_validate/RNA_temp2/recal/
bcftools query -f '%CHROM\t%POS\t%INFO\n' recal.simplified_${chr}.bcf > ${chr}_all.txt
}
main
echo "Job completed at $(date '+%d_%m_%y_%H_%M_%S')"


#=========================================================================#
#         Cleanup  DO NOT REMOVE OR CHANGE                                #
#=========================================================================#
cd $SLURM_SUBMIT_DIR
rm -rf /scratch/$USER/$SLURM_JOBID
