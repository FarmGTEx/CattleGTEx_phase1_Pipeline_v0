ARGS <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(data.table)
library(ggplot2)
####filter_samples and covert to bfile


vcf_file <- ARGS[1]	 #VCF file of filtered RNA-seq samples
reference_file <- ARGS[2]  #metadata containing the breed information
count_threshold <- ARGS[3]  #Minimum sample count threshold, 15 was used in our project
INFO_threshold <- ARGS[4]  #Minimum INFO value threshold, 0.95 was used in our project


INFO_file <- " "  #INFO value of each SNP
ref_file <- " "  #Merged reference panel
ref_meta <- " " #metadata of WGS samples containing the breed information
 
#extract sample names from VCF files
system(paste0("bcftools query -l ",vcf_file," |sed \"s/_bqsr//g;s/.bqsr//g\" > sample_list.txt")) 
system(paste0("bcftools query -l ",ref_file," > ref_list.txt"))
vcf_sample <- fread("sample_list.txt",header = F, data.table = F)[,1]
reference <- fread(reference_file,header = T ,data.table = F)

#match the breed info to sample names
ref_info <- fread(ref_meta) 
ref_sample <- fread("ref_list.txt",header = F, data.table = F)[,1]
ref_breed <- ref_info[ref_info$ID %in% ref_sample,c(1,3)]
select_reference <- reference[reference$BioSample %in% vcf_sample,]
data_breed <- as.data.frame(cbind(select_reference$BioSample,select_reference$Breed_Class))
colnames(data_breed) <- colnames(ref_breed)
all_data <- as.data.frame(rbind(ref_breed,data_breed))

#filter crossed species and breeds with sample less than 15 
summary <- as.data.frame(table(all_data$breed_new))
summary <- summary[order(summary$Freq,decreasing = T),]
summary <- summary[summary[,1]!= "Unknown",]
summary <- subset(summary, 
                  !grepl("Cross", summary[[1]]) & !grepl("x", summary[[1]]))
selected <- summary[summary$Freq >= count_threshold ,]
ref_extract <- ref_breed[ref_breed$breed_new %in% selected[,1],]
data_extract <- data_breed[data_breed$breed_new %in% selected[,1],]

write.table(ref_extract[,1],"ref_extract.txt",col.names = F , row.names = F ,quote = F)
write.table(unique(data_extract[,1]),"data_extract.txt",col.names = F , row.names = F ,quote = F)
sample_group <- rbind(ref_extract,unique(data_extract))
write.table(sample_group[,c(1,1,2)],"sample_group.txt",col.names = F , row.names = F ,quote = F,sep = "\t")

#match the INFO value to selected SNPs
system(paste0("cat ",INFO_file," | sed 's|=|\t|g ; s|;|\t|g' | cut -f 1,2,4,6 > INFO_file"))
INFO <- fread("INFO_file",header = F, data.table = F)
colnames(INFO) <- c("chr","pos","AF","INFO")
selected_SNP <- INFO[INFO$INFO >= INFO_threshold & INFO$AF > 0.05 & INFO$AF < 0.95, ]
write.table(selected_SNP[,1:2],paste0("SNP_extract_",INFO_threshold,".txt"),col.names = F , row.names = F ,quote = F,sep ="\t")

#extract samples and SNPs from the WGS reference panel and RNA-imputed VCF, and convert them into PLINK binary format (bfile)
system(paste0("bcftools view -R SNP_extract_",INFO_threshold,".txt -S ref_extract.txt -Oz -o ref_select.vcf.gz /faststorage/project/farmgtex/gtex/cattle_ref/CattleSNPs/panel.vcf.gz"))
system(paste0("bcftools view -R SNP_extract_",INFO_threshold,".txt -S data_extract.txt -Oz -o data_select.vcf.gz ",vcf_file))
system("bcftools index -f ref_select.vcf.gz")
system("bcftools index -f data_select.vcf.gz")
system("bcftools merge data_select.vcf.gz ref_select.vcf.gz -Ov -o impute_samples.vcf")
system("plink --vcf impute_samples.vcf --make-bed --out impute_samples --cow --double-id")
system("awk '{ $2 = $1 \"_\" $4; print }' impute_samples.bim > temp.bim && mv temp.bim impute_samples.bim")






#Calculate AED and fst
system("sed -i 's/ /_/g' sample_group.txt")
system("gcta64 --bfile impute_samples --fst --autosome-num 29 --sub-popu sample_group.txt --out fst_test")
result <- fread("fst_test.fst",header = T, data.table = F)
info_df <- result[,5:(dim(result)[2]-1)]
rownames(info_df)<-result[,2]
freq_matrix <- as.matrix(info_df)
average_euclidean_distances <- apply(freq_matrix, 1, function(row) {
  dist_matrix <- dist(row, method = "euclidean")
  mean(dist_matrix)
})
result$AED <- average_euclidean_distances
data_save <- result[,(dim(result)[2]-1):dim(result)[2]]
rownames(data_save)<-result[,2]
write.table(data_save,"FST_AED.txt",col.names = F , row.names = T ,quote = F,sep ="\t")
