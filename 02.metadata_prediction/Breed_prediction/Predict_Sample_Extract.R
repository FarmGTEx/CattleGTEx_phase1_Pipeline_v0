ARGS <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(data.table)
library(ggplot2)
####filter_samples and covert to bfile

vcf_file <- ARGS[1]	 #VCF file of RNA-seq samples
reference_file <- ARGS[2]  #metadata containing the breed information
count_threshold <- ARGS[3]  #Minimum sample count threshold, 15 was used in our project
INFO_threshold <- ARGS[4]  #Minimum INFO value threshold, 0.95 was used in our project

INFO_file <- " "  #INFO value of each SNP
ref_file <- " "  #Merged reference panel
ref_meta <- " " #metadata of WGS samples containing the breed information

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
summary <- as.data.frame(table(all_data$breed_new))
unknown_count <- summary[summary$Var1 == "Unknown", "Freq"]

other_count <- sum(summary$Freq[summary$Var1 != "Unknown"])

cat("Unknown count: ", unknown_count, "\n")
cat("Other breeds count: ", other_count, "\n")



summary <- summary[order(summary$Freq,decreasing = T),]
summary <- subset(summary, !grepl("Cross", summary[[1]]))
summary <- summary[summary[,1]!= "Unknown",]
summary <- subset(summary, 
                  !grepl("Cross", summary[[1]], ignore.case = TRUE) &
                  !grepl("x", summary[[1]], ignore.case = TRUE) &
                  !grepl("×", summary[[1]], ignore.case = TRUE))

Unknown_data <- data_breed[data_breed$breed_new =="Unknown",]
colnames(Unknown_data) <- colnames(ref_breed)
selected <- summary[summary$Freq >= count_threshold ,] 


ref_extract <- ref_breed[ref_breed$breed_new %in% selected[,1],]
data_extract <- data_breed[data_breed$breed_new %in% selected[,1],]
data_extract <- rbind(data_extract,Unknown_data)
write.table(unique(data_final),"predict_data_extract.txt",col.names = F , row.names = F ,quote = F)


system(paste0("bcftools view -R SNP_extract_",INFO_threshold,".txt -S predict_data_extract.txt -Oz -o predict_data_select.vcf.gz ",vcf_file))
system("bcftools index predict_data_select.vcf.gz")
system("bcftools merge predict_data_select.vcf.gz ref_select.vcf.gz -Ov -o predict_samples.vcf")
system("plink --vcf predict_samples.vcf --make-bed --out predict_samples --cow --double-id")
system("awk '{ $2 = $1 \"_\" $4; print }' predict_samples.bim > temp.bim && mv temp.bim predict_samples.bim")