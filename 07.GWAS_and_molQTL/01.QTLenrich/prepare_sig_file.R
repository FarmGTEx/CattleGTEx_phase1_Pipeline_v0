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

##prepare inout files for QTLenrich
path <- "/faststorage/project/farmgtex/GWAS_coloc/rank_QTL/isoQTL/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0("/faststorage/project/farmgtex/QTL_result_new/isoQTL/", i, "/")
}

df <- read.table("/faststorage/project/farmgtex/pipeline/isoQTL/gene_table", sep = "\t")
names(df) <- c("gene_id","pheno_id")

for (i in 11:length(tissue)) {
  setwd(paste0(path[[i]]))
  qtl <- fread(paste0(tissue[[i]], ".cis_qtl.txt.sig.gz"))
  #qtl <- qtl[qtl$qval_g1 < 0.05 & qtl$pval_g1 < qtl$pval_g1_threshold,]
  all_nom_qtl <- data.frame()
  for (chr in 1:29) {
    chr_qtl <- fread(paste0(tissue[[i]], ".cis_qtl_pairs.", chr, ".txt.gz"))
    #chr_qtl <- chr_qtl[chr_qtl$pval_g1 < 0.05, ]
    all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
  }
  all_nom_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% qtl$pheno_id,]
  qtl <- qtl[,c(3,14)]
  all_nom_qtl <- all_nom_qtl %>% left_join(qtl, by = "pheno_id")
  all_nom_qtl <- all_nom_qtl[all_nom_qtl$pval_g1 < all_nom_qtl$pval_g1_threshold,]
  qtl <- all_nom_qtl[,c(2,1)]
  names(qtl)[1] <- "variant_id"
  qtl$variant_id <- paste0("chr", qtl$variant_id)
  #qtl$gene_id <- str_extract(qtl$pheno_id, "ENSBTAG[0-9]+")
  qtl <- qtl %>% inner_join(df, by = "pheno_id")
  qtl <- qtl[,c(1,3)]
  #all_qtl <- rbind(all_qtl, qtl)
  #all_qtl <- all_qtl[!duplicated(all_qtl[, c("variant_id", "gene_id")]), ]
  fwrite(qtl, paste0("/faststorage/project/farmgtex/GWAS_coloc/QTLenrich/isoQTL/sig_pairs/", tissue[[i]], ".cis_qtl_pairs.significant.txt"), sep = "\t")
}

##prepare inout files for QTLenrich
path <- "/faststorage/project/farmgtex/GWAS_coloc/rank_QTL/3aQTL/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0("/faststorage/project/farmgtex/QTL_result_new/3aQTL/", i, "/")
}

for (i in 1:length(tissue)) {
  setwd(paste0(path[[i]]))
  qtl <- fread(paste0(tissue[[i]], ".cis_qtl.txt.sig.gz"))
  #qtl <- qtl[qtl$qval_g1 < 0.05 & qtl$pval_g1 < qtl$pval_g1_threshold,]
  all_nom_qtl <- data.frame()
  for (chr in 1:29) {
    chr_qtl <- fread(paste0(tissue[[i]], ".cis_qtl_pairs.", chr, ".txt.gz"))
    #chr_qtl <- chr_qtl[chr_qtl$pval_g1 < 0.05, ]
    all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
  }
  all_nom_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% qtl$pheno_id,]
  qtl <- qtl[,c(1,12)]
  all_nom_qtl <- all_nom_qtl %>% left_join(qtl, by = "pheno_id")
  all_nom_qtl <- all_nom_qtl[all_nom_qtl$pval_g1 < all_nom_qtl$pval_g1_threshold,]
  qtl <- all_nom_qtl[,c(2,1)]
  names(qtl)[1] <- "variant_id"
  qtl$variant_id <- paste0("chr", qtl$variant_id)
  qtl$gene_id <- str_extract(qtl$pheno_id, "ENSBTAG[0-9]+")
  #qtl <- qtl %>% inner_join(df, by = "pheno_id")
  qtl <- qtl[,c(1,3)]
  #all_qtl <- rbind(all_qtl, qtl)
  #all_qtl <- all_qtl[!duplicated(all_qtl[, c("variant_id", "gene_id")]), ]
  fwrite(qtl, paste0("/faststorage/project/farmgtex/GWAS_coloc/QTLenrich/3aQTL/sig_pairs/", tissue[[i]], ".cis_qtl_pairs.significant.txt"), sep = "\t")
}

##prepare inout files for QTLenrich
path <- "/faststorage/project/farmgtex/GWAS_coloc/rank_QTL/sQTL/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0("/faststorage/project/farmgtex/QTL_result_new/sQTL/", i, "/")
}

df <- read.table("/faststorage/project/farmgtex/pipeline/sQTL/cluster_groups.txt", sep = "\t")
names(df) <- c("pheno_id","gene_id")

for (i in 1:length(tissue)) {
  setwd(paste0(path[[i]]))
  qtl <- fread(paste0(tissue[[i]], ".cis_qtl.txt.sig.gz"))
  #qtl <- qtl[qtl$qval_g1 < 0.05 & qtl$pval_g1 < qtl$pval_g1_threshold,]
  all_nom_qtl <- data.frame()
  for (chr in 1:29) {
    chr_qtl <- fread(paste0(tissue[[i]], ".cis_qtl_pairs.", chr, ".txt.gz"))
    #chr_qtl <- chr_qtl[chr_qtl$pval_g1 < 0.05, ]
    all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
  }
  all_nom_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% qtl$pheno_id,]
  qtl <- qtl[,c(3,14)]
  all_nom_qtl <- all_nom_qtl %>% left_join(qtl, by = "pheno_id")
  all_nom_qtl <- all_nom_qtl[all_nom_qtl$pval_g1 < all_nom_qtl$pval_g1_threshold,]
  qtl <- all_nom_qtl[,c(2,1)]
  names(qtl)[1] <- "variant_id"
  qtl$variant_id <- paste0("chr", qtl$variant_id)
  #qtl$gene_id <- str_extract(qtl$pheno_id, "ENSBTAG[0-9]+")
  qtl <- qtl %>% inner_join(df, by = "pheno_id")
  qtl <- qtl[,c(1,3)]
  #all_qtl <- rbind(all_qtl, qtl)
  #all_qtl <- all_qtl[!duplicated(all_qtl[, c("variant_id", "gene_id")]), ]
  fwrite(qtl, paste0("/faststorage/project/farmgtex/GWAS_coloc/QTLenrich/sQTL/sig_pairs/", tissue[[i]], ".cis_qtl_pairs.significant.txt"), sep = "\t")
}

##prepare inout files for QTLenrich
path <- "/faststorage/project/farmgtex/GWAS_coloc/rank_QTL/stQTL/"
tissue <- list.dirs(path, full.names = TRUE, recursive = FALSE)
tissue <- sapply(tissue, function(x) unlist(strsplit(x, "\\/"))[9])
tissue <- data.frame(tissue)
tissue <- tissue$tissue
path<-NULL
for(i in tissue){
  path[[i]]<-paste0("/faststorage/project/farmgtex/QTL_result_new/stQTL/", i, "/")
}

for (i in 1:length(tissue)) {
  setwd(paste0(path[[i]]))
  qtl <- fread(paste0(tissue[[i]], ".cis_qtl.txt.sig.gz"))
  #qtl <- qtl[qtl$qval_g1 < 0.05 & qtl$pval_g1 < qtl$pval_g1_threshold,]
  all_nom_qtl <- data.frame()
  for (chr in 1:29) {
    chr_qtl <- fread(paste0(tissue[[i]], ".cis_qtl_pairs.", chr, ".txt.gz"))
    #chr_qtl <- chr_qtl[chr_qtl$pval_g1 < 0.05, ]
    all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
  }
  all_nom_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% qtl$pheno_id,]
  qtl <- qtl[,c(1,12)]
  all_nom_qtl <- all_nom_qtl %>% left_join(qtl, by = "pheno_id")
  all_nom_qtl <- all_nom_qtl[all_nom_qtl$pval_g1 < all_nom_qtl$pval_g1_threshold,]
  qtl <- all_nom_qtl[,c(2,1)]
  names(qtl)[1] <- "variant_id"
  names(qtl)[2] <- "gene_id"
  qtl$variant_id <- paste0("chr", qtl$variant_id)
  #qtl$gene_id <- str_extract(qtl$pheno_id, "ENSBTAG[0-9]+")
  #qtl <- qtl %>% inner_join(df, by = "pheno_id")
  #qtl <- qtl[,c(1,3)]
  fwrite(qtl, paste0("/faststorage/project/farmgtex/GWAS_coloc/QTLenrich/stQTL/sig_pairs/", tissue[[i]], ".cis_qtl_pairs.significant.txt"), sep = "\t")
}