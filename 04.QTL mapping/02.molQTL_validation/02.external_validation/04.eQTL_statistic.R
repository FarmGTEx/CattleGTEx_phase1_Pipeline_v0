library(data.table)
library(dplyr)

source <- "RNA"
type <- "result"



result_path <- paste0("/faststorage/project/farmgtex/pipeline/QTL_validation2/",source,"/",type,"/")

prefix <- source



summary <- fread(paste0(result_path,prefix,".cis_qtl.txt.gz"),header = T, data.table = F)

egene <- summary[summary$pval_g1 < summary$pval_g1_threshold & summary$qval_g1 < 0.05,1]

final_result <- NULL
for (i in 1:29){
    temp_result <- fread(paste0(result_path,prefix,".cis_qtl_pairs.",i,".txt.gz"),header = T, data.table = F)
    temp_result <- temp_result[temp_result[,1] %in% egene,]
    temp_result$pval_g1_threshold <- summary$pval_g1_threshold[match(temp_result[,1],summary[,1])]
    result_with_fdr <- temp_result %>%
    group_by(pheno_id) %>%  # 按 pheno_id 分组
    mutate(fdr_pval_g1 = p.adjust(pval_g1, method = "fdr")) %>%  # FDR 校正
    ungroup()%>%
    as.data.frame()
    significant_snp <- result_with_fdr[result_with_fdr$pval_g1 < result_with_fdr$pval_g1_threshold & result_with_fdr$fdr_pval_g1 < 0.05,]
        if (i == 1) {
        final_result <- significant_snp  # 如果是第一个文件，直接赋值
    } else {
        final_result <- rbind(final_result, significant_snp)  # 否则追加
    }    
}
write.table(final_result[,1:7],paste0(result_path,prefix,".cis_qtl.txt.sig_all.gz"),quote = F,col.names = T,row.names= F,sep ="\t")
