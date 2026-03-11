# Calculate gene expression principal components (PCs) to use as covariates for correcting confounding factors within each tissue


args <- commandArgs(trailingOnly = TRUE)
type <- args[1]


library(data.table)
library(dendextend)
library(limma)

syscolor <- fread(" ") #color for system
tiscolor <- fread(" ") #color for tissue
germs <- fread(" ") #tissue-germ map
gemcolor <- fread(" ") #color for germs

reference <- fread(" ") #metadata of all RNA-seq samples
samples <- fread(" ") #sample list used in molQTL mapping


##match colors for each sample
matrix <- readRDS(paste0(type, "_cluster.rds"))
sample_names <- colnames(matrix)
ref <- reference[match(sample_names, reference$BioSample), ]

ref$germ <- germs[match(ref$Predict_Tissue, germs[,1]),2]
ref$Germ_color <- gemcolor[match(ref$germ, gemcolor[,1]), 2]
ref$Tissue_color <- tiscolor[match(ref$Predict_Tissue, tiscolor[,1]), 2]
ref$Sys_color <- syscolor[match(ref$Predict_System, syscolor[,1]), 2]
matrix <- matrix[, ref$BioSample, drop = FALSE]

pcs <- readRDS("PCs.rds")
expr_log <- log2(matrix + 0.01)

##Correct molecular phenotype matrix using limma with gene expression PCs as covariates
correct_within_tissue <- function(expr, tissue_df, pcs) {
  
  corrected <- expr
  tissues <- unique(tissue_df$Predict_Tissue)
  
  for (tissue in tissues) {
    
    tissue_samples <- tissue_df$BioSample[tissue_df$Predict_Tissue == tissue]
    common_samples <- intersect(tissue_samples, colnames(expr))
    
    if (length(common_samples) >= 40) {
      
      tissue_expr <- expr[, common_samples, drop = FALSE]
      
      gene_vars <- apply(tissue_expr, 1, var, na.rm = TRUE)
      valid_genes <- which(gene_vars > 0 & !is.na(gene_vars))
      
      if (length(valid_genes) == 0) next
      
      tissue_expr_filtered <- tissue_expr[valid_genes, , drop = FALSE]
      
      tryCatch(
        {
          ## PCA covariates
          pca_cov <- as.matrix(pcs[[tissue]])
          pca_cov <- pca_cov[common_samples, , drop = FALSE]
          
          tissue_corrected_filtered <- removeBatchEffect(
            tissue_expr_filtered,
            covariates = pca_cov
          )
          
          corrected[valid_genes, common_samples] <- tissue_corrected_filtered
        },
        error = function(e) {
          message(
            "Error in tissue ", tissue, ": ", e$message
          )
        }
      )
    }
  }
  
  return(corrected)
}



expr_corrected <- correct_within_tissue(expr_log, ref, pcs)

data_corrected <- t(expr_corrected)
data_scaled <- scale(data_corrected)
dist_corrected <- dist(data_scaled, method = "manhattan")
hc_corrected <- hclust(dist_corrected)
dend_corrected <- as.dendrogram(hc_corrected)

saveRDS(dist_corrected, paste0(type, "_corrected_dist.rds"))
saveRDS(dend_corrected, paste0(type, "_corrected_dend.rds"))


### plot the corrected clustering
pdf(paste0(type, "_corrected_cluster.pdf"), width = 10, height = 6)
par(mar = c(10, 6, 6, 8))

tree_h <- attr(dend_corrected, "height")
bar_h  <- tree_h * bar_ratio

plot(dend_corrected, type = "rectangle", leaflab = "none", axes = FALSE,
     edgePar = list(col = "gray40", lty = 1, lwd = 2))

colored_bars(colors = ref$Tissue_color, dend = dend_corrected, 
             y_shift = -bar_h * 0.2, # Shift slightly below tree
             y_scale = bar_h, 
             rowLabels = "Tissue")

colored_bars(colors = ref$Sys_color, dend = dend_corrected,
             y_shift = -(bar_h * 1.7), # Shift below the first bar
             y_scale = bar_h, 
             rowLabels = "System",
             add = TRUE)
             
colored_bars(colors = ref$Germ_color, dend = dend_corrected,
             y_shift = -(bar_h * 3.2), # Shift below the first bar
             y_scale = bar_h, 
             rowLabels = "Germ",
             add = TRUE)             

dev.off()


data_original <- t(expr_log)
data_original_scaled <- scale(data_original)
dist_original <- dist(data_original_scaled, method = "manhattan")
hc_original <- hclust(dist_original)
dend_original <- as.dendrogram(hc_original)

saveRDS(dist_original, paste0(type, "_original_dist.rds"))
saveRDS(dend_original, paste0(type, "_original_dend.rds"))


### plot the original clustering
pdf(paste0(type, "_original_cluster.pdf"), width = 10, height = 6)
par(mar = c(10, 6, 6, 8))

# Get tree height for this specific dendrogram
tree_h_orig <- attr(dend_original, "height")
bar_h_orig  <- tree_h_orig * bar_ratio 

plot(dend_original, type = "rectangle", leaflab = "none", axes = FALSE,
     edgePar = list(col = "gray40", lty = 1, lwd = 2))

colored_bars(colors = ref$Tissue_color, dend = dend_original, 
             y_shift = -bar_h_orig * 0.2, 
             y_scale = bar_h_orig, 
             rowLabels = "Tissue")

colored_bars(colors = ref$Sys_color, dend = dend_original,
             y_shift = -(bar_h_orig * 1.7), 
             y_scale = bar_h_orig, 
             rowLabels = "System",
             add = TRUE)
             
colored_bars(colors = ref$Germ_color, dend = dend_original,
             y_shift = -(bar_h_orig * 3.2), # Shift below the first bar
             y_scale = bar_h_orig, 
             rowLabels = "Germ",
             add = TRUE)  
dev.off()

