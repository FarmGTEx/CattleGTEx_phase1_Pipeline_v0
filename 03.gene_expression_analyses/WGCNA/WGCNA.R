# Apply WGCNA to raw gene expression counts for each tissue

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
library(PCAForQTL)
library(WGCNA)
library(data.table)
library(dendextend)
library(limma)
library(edgeR)
#---------------------------------
 Convert counts to TMM-normalized values and correct using expression PCs
#---------------------------------
count <- readRDS(tissue,"/count_QTL.rds"))
dge <- DGEList(counts = count)
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

gene_vars <- apply(logCPM, 1, var, na.rm = TRUE)
valid_genes <- which(gene_vars > 0 & !is.na(gene_vars))
tissue_expr_filtered <- logCPM[valid_genes, , drop = FALSE]

pca <- prcomp(t(tissue_expr_filtered), center = TRUE, scale. = TRUE)
PCs <- pca$x

elbow <- PCAForQTL::runElbow(prcompResult = pca)
pca_results <- PCs[,1:elbow]

tissue_corrected_filtered <- removeBatchEffect(
    tissue_expr_filtered,
    covariates = pca_results
)

#---------------------------------
 Perform WGCNA and save the results
#---------------------------------
datExpr <- tissue_corrected_filtered
gene_vars <- apply(datExpr, 1, var, na.rm = TRUE)
valid_genes <- which(gene_vars > 0 & !is.na(gene_vars))
datExpr <- datExpr[valid_genes, , drop = FALSE]

datExpr <- t(datExpr)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()  # ���ö��߳�

net <- blockwiseModules(
  datExpr,
  TOMType = "unsigned",       # �޷�������
  reassignThreshold = 0,
  mergeCutHeight = 0.25,      # �ϲ�ģ����ֵ
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  corType = "pearson",
  verbose = 3
)

resfile <- paste0("WGCNA_", tissue, ".RData")
saveRDS(net, file = file.path('/faststorage/project/farmgtex/pipeline/WGCNA/1_RData', resfile))

####
moduleColors <- labels2colors(net$colors)
Sta <- data.frame(table(moduleColors))
Sta$Module.num <- length(Sta$moduleColors)
write.csv(Sta,
          file.path('/faststorage/project/farmgtex/pipeline/WGCNA/2_separate_gene.num',
                    paste0("WGCNA_", tissue, "_module_gene_count.csv")),
          quote = FALSE, row.names = FALSE)        
          
####
res_mat <- data.frame(net$MEs)
sample_order <- rownames(datExpr)
rownames(res_mat) <- sample_order
res_mat <- res_mat[, !colnames(res_mat) %in% "MEgrey"]  

write.csv(res_mat,
          file.path('/faststorage/project/farmgtex/pipeline/WGCNA/3_separate_eigen',
                    paste0("WGCNA_", tissue, "_eigengenes.csv")),
          quote = FALSE)
          
####
gene_name <- data.frame(gene_id = colnames(datExpr))
moduleLabels <- data.frame(net.colors = paste0("ME", net$colors))
moduleColors_df <- data.frame(labels2colors(net$colors))

geneInfo0 <- cbind.data.frame(
  "gene" = colnames(datExpr),
  gene_name,
  moduleLabels,
  moduleColors_df
)

names(geneInfo0) <- c("gene", "gene_id", "module_MEnumber", "module_Color")
geneInfo0 <- geneInfo0[geneInfo0$module_Color != "grey", ]

write.table(geneInfo0,
            file.path('/faststorage/project/farmgtex/pipeline/WGCNA/4_separate_clusters',
                      paste0("WGCNA_", tissue, "_geneInfo.tsv")),
            quote = FALSE, sep = "\t", row.names = FALSE)

MEs <- net$MEs
MEs <- MEs[, !colnames(MEs) %in% "MEgrey"]  
kME_all <- signedKME(datExpr, MEs)  
geneInfo0$module_num <- as.numeric(sub("ME", "", geneInfo0$module_MEnumber))
colnames(kME_all) <- as.numeric(sub("kME", "", colnames(kME_all)))

gene_kME <- data.frame(
  gene = geneInfo0$gene,
  gene_id = geneInfo0$gene_id,
  module_Color = geneInfo0$module_Color,
  module_MEnumber = geneInfo0$module_MEnumber,
  kME = mapply(function(g, mod){
    kME_all[match(g, rownames(kME_all)), match(mod, colnames(kME_all))]
  }, geneInfo0$gene, geneInfo0$module_num)
)

outdir_kME <- '/faststorage/project/farmgtex/pipeline/WGCNA/5_gene_kME'
write.csv(gene_kME,
          file.path(outdir_kME, paste0("WGCNA_", tissue, "_gene_kME.csv")),
          quote = FALSE, row.names = FALSE)
     