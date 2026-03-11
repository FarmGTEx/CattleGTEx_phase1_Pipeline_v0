options(stringsAsFactors = FALSE)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
library(dplyr)
library(coloc)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)

path <- "/faststorage/project/farmgtex/GWAS_coloc/rank_QTL/eQTL/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0("/faststorage/project/farmgtex/QTL_result_new/eQTL/", i)
}
gtf = rtracklayer::import("/faststorage/project/cattle_gtexs/reference/Bos_taurus.ARS-UCD1.2.110.gtf")
gtf = as.data.frame(gtf)
gtf <- gtf[gtf$type == "gene",]
gtf <- gtf[,c(5,10)]
names(gtf) <- c("Orientation","pheno_id")

for (i in 16:length(tissue)) {
  setwd(path[[i]])
  all_nom_qtl <- data.frame()
  for (chr in 1:29) {
    chr_qtl <- fread(paste0(tissue[[i]], ".cis_qtl_pairs.", chr, ".txt.gz"))
    #chr_qtl <- chr_qtl[chr_qtl$pval_g1 < 0.05, ]
    all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
  }
  all_nom_qtl$t_stat <- all_nom_qtl$beta_g1/all_nom_qtl$beta_se_g1
  all_nom_qtl$FDR <- p.adjust(all_nom_qtl$pval_g1, method = "BH")
  all_nom_qtl <- all_nom_qtl[,c(2,1,5,8,7,9)]
  names(all_nom_qtl) <- c("SNP","gene","beta","t-stat","p-value","FDR")
  #all_nom_qtl$Chr <- sub("_.*", "", all_nom_qtl$variant_id)
  #all_nom_qtl$BP <- sub("^[^_]*_([^_]*)_.*", "\\1", all_nom_qtl$variant_id)
  #all_nom_qtl$A1 <- sub("^[^_]*_[^_]*_([^_]*)_.*", "\\1", all_nom_qtl$variant_id)
  #all_nom_qtl$A2 <- sub(".*_([^_]*)$", "\\1", all_nom_qtl$variant_id)
  #all_nom_qtl <- all_nom_qtl %>% inner_join(gtf, by = "pheno_id")
  #all_nom_qtl <- all_nom_qtl[,c(2,8,9,10,11,4,1,8,9,1,12,5,6,7)]
  #all_nom_qtl$variant_id <- paste0(all_nom_qtl$Chr, "_", all_nom_qtl$BP)
  #names(all_nom_qtl) <- c("SNP","Chr","BP","A1","A2","Freq","Probe","Probe_Chr","Probe_bp","Gene","Orientation","b","se","p")
  write.table(all_nom_qtl, paste0("/faststorage/project/farmgtex/GWAS_coloc/SMR/all_pairs/", tissue[[i]], ".eqtl_allpairs.txt"), sep = "\t", quote = F, row.names = F)
}