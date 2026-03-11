find /faststorage/project/farmgtex/fastq_final/ -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sed 's/$/_bqsr.bam/'  > all_sample
ls "/faststorage/project/farmgtex/gtex/bqsr_bamlist/" | grep -f all_sample | grep -v "SA" > /faststorage/project/farmgtex/gtex/bqsr_bamlist/RNA_list

####
cat /faststorage/project/farmgtex/gtex/cattle_ref/SNPs/split/bins.txt|while read file
do
sbatch RNA_geno.sh $file 1 3
done

####
cat /faststorage/project/farmgtex/gtex/cattle_ref/chrlist |while read file
do
sbatch RNA_geno.sh 1 $file 4
done

###
sbatch RNA_geno.sh 1 1 5
###
sbatch RNA_geno.sh 1 1 7


#####
ls /faststorage/project/farmgtex/gtex/bqsr_WGS/*vcf.gz | xargs -n 1 basename | sed "s/\.g\.vcf\.gz//g" > WGS_sample
cat all_sample | sed "s/_bqsr\.bam//g" > RNA_sample
cat WGS_sample | grep -f RNA_sample  | sed 's/$/.g.vcf.gz/' | awk '{print "--variant "$0}' > gvcf_args.txt


cat WGS_sample | grep -f RNA_sample  | sed 's/$/.g.vcf.gz/' | while read sample
do
gatk IndexFeatureFile -I /faststorage/project/farmgtex/gtex/bqsr_WGS/${sample}
done