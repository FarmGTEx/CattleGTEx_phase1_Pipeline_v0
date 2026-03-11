library(data.table)
ARGS <- commandArgs(trailingOnly = TRUE)
scenario_select = ARGS[1]
dapars2_result_file <- paste0("/faststorage/project/farmgtex/pipeline/polyA/",scenario_select,"/",scenario_select,"_final_predict_clean.txt")
gene_gtf_file <- "/faststorage/project/farmgtex/pipeline/eQTL/TSS.gtf"
gtf <- fread(gene_gtf_file,header = T,data.table = F)
result <- read.table(dapars2_result_file,header = T)
split_genes <- strsplit(result$Gene, "\\|")
gene_columns <- as.data.frame(do.call(rbind, split_genes))
final_df <- cbind(gene_columns, result[,1]) 
colnames(final_df) <- c("transcript_id", "gene_id", "chromosome", "strand","PolyAloci")
merged_df <- merge(final_df, gtf, by = "gene_id")
write.table(merged_df[,c(6:8,5)],paste0("/faststorage/project/farmgtex/pipeline/3aQTL/",scenario_select,"/annotation.gtf"),quote = F, sep = "\t" , col.names= T, row.names = F)