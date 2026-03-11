# Calculate gene expression principal components (PCs) to use as covariates for correcting confounding factors within each tissue

library(data.table)
library(dendextend)
library(limma)
library(PCAForQTL)

type="gene"

reference <- fread(" ") #metadata for RNA-seq samples
samples <- fread(" ") #sample list used in molQTL mapping

#
matrix <- readRDS(paste0(type, "_cluster.rds"))
sample_names <- colnames(matrix)
ref <- reference[match(sample_names, reference$BioSample), ]
matrix <- matrix[, ref$BioSample, drop = FALSE]

pca_results <- list()

for (tis in keep_tissue) {

  message("Processing tissue: ", tis)

  idx <- which(ref$Predict_Tissue == tis)
  samples_tis <- ref$BioSample[idx]

  if (length(samples_tis) < 10) {
    message("  Skip (n < 10)")
    next
  }

  mat_tis <- as.data.frame(t(matrix[, samples_tis, drop = FALSE]))

  mat_tis2 <- mat_tis[, apply(mat_tis, 2, sd) != 0]

  pca <- prcomp(mat_tis2, center = TRUE, scale. = TRUE)
  PCs<-pca$x
  elbow <- PCAForQTL::runElbow(
    prcompResult= pca
  )
  pca_results[[tis]] <- PCs[,1:elbow]
}

saveRDS(pca_results,"PCs.rds")
