library(data.table)
ARGS <- commandArgs(trailingOnly = TRUE)
QTL_scenarios <- fread("/faststorage/project/farmgtex/pipeline/QTL_design",header = F, data.table =F)
scenario_select = ARGS[1]
scenario <- gsub("_"," ", scenario_select)
tis <- QTL_scenarios[QTL_scenarios[,1] == scenario,2]

sample <- fread("/faststorage/project/farmgtex/gtex/VCFfile/chr_final/recal/final_vcf/deduplicate/IBS_grm/uni_sample_list.txt",header = T, data.table =F)
sample <- sample[sample$Tissue %in% tis,]
leave <- unique(sample$leave)
count <- readRDS("/faststorage/project/farmgtex/pipeline/enhancer/enhancer_count_matrix_new.rds")

select_count <- count[,colnames(count) %in% leave]
saveRDS(select_count,paste0("/faststorage/project/farmgtex/pipeline/enQTL/",scenario_select,"/count_QTL.rds"))

write.table(colnames(select_count),paste0("/faststorage/project/farmgtex/pipeline/enQTL/",scenario_select,"/sample_list.txt"),quote = F, sep = "\t" , col.names= F, row.names = F)