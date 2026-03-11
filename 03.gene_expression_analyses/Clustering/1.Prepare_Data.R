library(data.table)

reference <- fread(" ") #metadata of all RNA-seq samples
samples <- fread(" ") #sample list used in molQTL mapping

tpm_threshold <- 0.1

###gene expression/ isoform expression/ exon expression/ enhancer expression/

# extract the top phenotypes with highest standard deviation
matrix <- readRDS(" ") ##the molecular phenotype expression matrix
sub_matrix <- matrix[,colnames(matrix)%in% samples]
gene_sd <- apply(sub_matrix, 1, sd, na.rm = TRUE)
top5000_genes <- names(sort(gene_sd, decreasing = TRUE))[1:5000]
top5000_matrix <- sub_matrix[top5000_genes, ]

saveRDS(top5000_matrix,"cluster.rds")

###alternative splicing
# Extract top intron clusters with highest standard deviation, and remove clusters with zero counts in all samples
matrix <- readRDS("") #psi matrix for all samples 
matrix <- matrix[, colnames(matrix) %in% samples]
matrix<- matrix[rowSums(matrix)>0,] 
matrix$sd <- apply(matrix, 1, sd)
matrix <- matrix[order(matrix$sd, decreasing = T),]
matrix <- head(matrix, n = 5000)
matrix<- matrix[,-ncol(matrix)]
saveRDS(matrix,"splicing_matrix.rds")

###3' alternative polyadenylation
# Extract top loci with highest standard deviation, keeping only loci with at least 10 non-missing values

expr <- fread(" ") #PDUI matrix for all samples 
rownames(expr) <- expr[[1]]
expr <- expr[, -c(1:4)]
sample_names <- colnames(expr)

sample_names <- basename(sample_names)
sample_names <- sub("\\.bam_PDUI$", "", sample_names)
colnames(expr) <- sample_names
expr <- expr[, colnames(expr) %in% samples]

n_non_na <- apply(expr, 1, function(x) sum(!is.na(x)))
expr$sd <- apply(expr, 1, sd, na.rm = TRUE)

expr <- expr[
  n_non_na >= 10 & expr$sd > 0 & !is.na(expr$sd),
]
matrix <- expr [order(expr$sd, decreasing = T),]
matrix <- head(matrix, n = 5000)
matrix<- matrix[,-ncol(matrix)]
saveRDS(matrix,"3a_matrix.rds")

###RNA stability
# Extract top 5000 genes with highest standard deviation of exon/intron ratio, keeping genes with ≥10 non-missing values 

expr <- fread(" ") #exon/intron constitutive ratio matrix for all samples 
expr <- expr[, colnames(expr) %in% samples]
n_non_na <- apply(expr, 1, function(x) sum(!is.na(x)))
expr$sd <- apply(expr, 1, sd, na.rm = TRUE)
expr <- expr[
  n_non_na >= 10 & expr$sd > 0 & !is.na(expr$sd),
]
matrix <- expr [order(expr$sd, decreasing = T),]
matrix <- head(matrix, n = 5000)
matrix<- matrix[,-ncol(matrix)]