## Perform GO term analysis of tissue-specific genes, considering only genes with logFC > 1.5

library(data.table)
library(org.Bt.eg.db)
library(clusterProfiler)
library(dplyr)
library(AnnotationHub)


tissues <- fread(" ",data.table= F,header = F)[,1] ##eQTL tissue list


logFC_cutoff<-1.5

all_GO <- list() 

for (tissue in tissues){
    DEG <-read.table(paste0(tissue,"_Wilcoxon.rst.tsv"),header = T,stringsAsFactors = F)
    DEG <- DEG[DEG[,1] > 1.5 & DEG[,3] < 0.05, ]
    colnames(DEG) <- c("Log2FC","P_val","FDR")
    upgene<- DEG
    upgene<-upgene[order(upgene$FDR,decreasing = F),]
    gene_list<-rownames(upgene)
    ego <- enrichGO(
      gene         = gene_list,
      OrgDb        = org.Bt.eg.db,
      keyType      = "ENSEMBL",    
      ont          = "BP",          # Biological Process
      pAdjustMethod= "BH",
      pvalueCutoff = 0.05,
      readable     = TRUE           
    )
    
    if(!is.null(ego) && nrow(as.data.frame(ego)) > 0){
    ego_df <- as.data.frame(ego)[,c("ID","Description","p.adjust")]
    colnames(ego_df) <- c("ID","Description","FDR")
    ego_df$tissue <- tissue
    all_GO[[tissue]] <- ego_df
} else {
    message(paste("No GO term found for tissue:", tissue))
}
}

  GO_long <- bind_rows(all_GO)
  write.csv(GO_long, "GO_all_tissues_long.csv", row.names = F)