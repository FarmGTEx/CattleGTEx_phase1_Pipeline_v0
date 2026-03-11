ARGS <- commandArgs(trailingOnly = TRUE)
scenario_select = ARGS[1]


library(data.table)
"%&%" = function(a, b) { paste0(a, b) }
QTL_scenarios <- fread("/faststorage/project/farmgtex/pipeline/QTL_design",header = F, data.table =F)
scenario <- gsub("_"," ", scenario_select)
tis <- QTL_scenarios[QTL_scenarios[,1] == scenario,2]
sample <- fread("/faststorage/project/farmgtex/gtex/VCFfile/chr_final/recal/final_vcf/deduplicate/IBS_grm/uni_sample_list.txt",header = T, data.table =F)
sample <- sample[sample$Tissue %in% tis,]
leave <- unique(sample$leave)

#count <- fread("/faststorage/project/farmgtex/pipeline/splicing/cattle_30_perind_numers.counts",header = F,data.table = F)
#IDs <- fread("/faststorage/project/farmgtex/pipeline/splicing/IDs",data.table = F, header = F)
#rownames(count) <- count[,1]
#count <- count[,-1]
#colnames(count) <- IDs[1,]
#saveRDS(count,"/faststorage/project/farmgtex/pipeline/sQTL/count.rds")
count <- readRDS("/faststorage/project/farmgtex/pipeline/sQTL/count.rds")
select_count <- count[,colnames(count) %in% leave]

#filter out clusters in less than 50% samples 

non_zero_counts <- rowSums(select_count != 0)
threshold1 <- ncol(select_count) * 0.5

cluster_counts <- rowSums(select_count)
threshold2 <- max(10, 0.1 * ncol(select_count))

filtered <- select_count[non_zero_counts > threshold1&cluster_counts > threshold2, ]

#filter out clusters with low complexity

z_scores <- t(scale(t(filtered))) 

low_variability <- rowSums(abs(z_scores) < 0.25)


high_outliers <- rowSums(abs(z_scores) > 6)

n <- ncol(filtered)

remove <- filtered[
  low_variability >= (n - 3) & high_outliers <= 3, 
]

final <- filtered[!rownames(filtered) %in% rownames(remove),]

write.table(colnames(final),paste0("/faststorage/project/farmgtex/pipeline/sQTL/",scenario_select,"/sample_list.txt"),quote = F, sep = "\t" , col.names= F, row.names = F)
write.table(rownames(final),paste0("/faststorage/project/farmgtex/pipeline/sQTL/",scenario_select,"/cluster_list.txt"),quote = F, sep = "\t" , col.names= F, row.names = F)
