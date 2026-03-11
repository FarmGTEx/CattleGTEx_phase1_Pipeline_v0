library(data.table)
ARGS <- commandArgs(trailingOnly = TRUE)
QTL_scenarios <- fread("/faststorage/project/farmgtex/pipeline/QTL_design",header = F, data.table =F)
scenario_select = ARGS[1]
scenario <- gsub("_"," ", scenario_select)
tis <- QTL_scenarios[QTL_scenarios[,1] == scenario,2]
set.seed(123)  # 为了可重复的随机结果
sample <- fread("/faststorage/project/farmgtex/gtex/VCFfile/chr_final/recal/final_vcf/deduplicate/IBS_grm/uni_sample_list.txt",header = T, data.table =F)
sample <- sample[sample$Tissue %in% tis,]
leave <- unique(sample$leave)
TPM <- readRDS("/faststorage/project/farmgtex/pipeline/gene_TPM/cattleTPMmatrix.rds")
count <- readRDS("/faststorage/project/farmgtex/pipeline/gene_count/cattle_gene_count_matrix.rds")
half_n <- floor(length(leave) / 2)
discovery <- leave[1:half_n]
validation <- leave[(half_n + 1):length(leave)]


select_TPM_discovery <- TPM[, colnames(TPM) %in% discovery]
saveRDS(select_TPM_discovery, paste0("/faststorage/project/farmgtex/pipeline/internal_validation/", scenario_select, "/discovery/TPM_QTL.rds"))
select_count_discovery <- count[, colnames(count) %in% discovery]
saveRDS(select_count_discovery, paste0("/faststorage/project/farmgtex/pipeline/internal_validation/", scenario_select, "/discovery/count_QTL.rds"))

select_TPM_validation <- TPM[, colnames(TPM) %in% validation]
saveRDS(select_TPM_validation, paste0("/faststorage/project/farmgtex/pipeline/internal_validation/", scenario_select, "/validation/TPM_QTL.rds"))
select_count_validation <- count[, colnames(count) %in% validation]
saveRDS(select_count_validation, paste0("/faststorage/project/farmgtex/pipeline/internal_validation/", scenario_select, "/validation/count_QTL.rds"))