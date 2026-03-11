library(data.table)
ARGS <- commandArgs(trailingOnly = TRUE)
gene_filter <- "/faststorage/project/farmgtex/pipeline/genechosen"
QTL_scenarios <- fread("/faststorage/project/farmgtex/pipeline/QTL_design",header = F, data.table =F)
filter = fread(gene_filter,header = F,data.table = F)
scenario_select = ARGS[1]
scenario <- gsub("_"," ", scenario_select)
tis <- QTL_scenarios[QTL_scenarios[,1] == scenario,2]
sample <- fread("/faststorage/project/farmgtex/gtex/VCFfile/chr_final/recal/final_vcf/deduplicate/IBS_grm/uni_sample_list.txt",header = T, data.table =F)
sample <- sample[sample$Tissue %in% tis,]
leave <- unique(sample$leave)
exon <- readRDS("/faststorage/project/farmgtex/pipeline/RNA_stability/RNA_stability_exon.rds")
intron <- readRDS("/faststorage/project/farmgtex/pipeline/RNA_stability/RNA_stability_intron.rds")
exon <- exon[,colnames(exon) %in% leave]
intron <- intron[,colnames(intron) %in% leave]
intron <- intron[rownames(exon), ]
intron <- intron[,colnames(exon)]
ratio_matrix <- exon / intron
ratio_matrix <- ratio_matrix[rownames(ratio_matrix) %in% filter[,1],]
ratio_matrix[is.infinite(as.matrix(ratio_matrix))] <- NA
ratio_matrix[is.nan(as.matrix(ratio_matrix))] <- NA
na_ratio <- rowSums(is.na(ratio_matrix)) / ncol(ratio_matrix)
ratio_matrix <- ratio_matrix[na_ratio <= 0.5, ]

write.table(ratio_matrix,paste0("/faststorage/project/farmgtex/pipeline/stQTL/",scenario_select,"/rna_ratio.txt"),quote = F, sep = "\t" , col.names=T, row.names = T)
write.table(colnames(ratio_matrix),paste0("/faststorage/project/farmgtex/pipeline/stQTL/",scenario_select,"/sample_list.txt"),quote = F, sep = "\t" , col.names= F, row.names = F)