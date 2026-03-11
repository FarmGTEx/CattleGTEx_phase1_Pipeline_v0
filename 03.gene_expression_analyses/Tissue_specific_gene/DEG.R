##Apply the Wilcoxon test to identify tissue-specific expressed genes

ARGS <- commandArgs(trailingOnly = TRUE)
tissue_name <- ARGS[1]
tissue <- gsub("_"," ",tissue_name)
TPM_file <- " " #TPM matrix of all samples 
reference_file <- " " #metadata of all samples
tissue_file <- " " #tissue list (higher than 40 samples)
filter_file <- " " #sample list

DEG_file <- " " #tissue group file
#e.g.
#Immune System   Tonsil
#Immune System   Blood
#Gut     Jejunum
#Gut     Ileum


library(data.table)
suppressWarnings(library(edgeR, quietly = T))


####read files 
matrix <- readRDS(TPM_file)
reference <- fread(reference_file,header = T, data.table = F)
tissue_info <- fread(tissue_file,header = F, data.table = F)
filter <- fread(filter_file,header = F, data.table = F)
group <- fread(DEG_file, header = F, data.table = F)

####extract samples (filter) and genes (expressed) 
sample <- intersect(filter[,1],reference[reference[,6] %in% tissue_info$V1,1])
tissue_matrix <- matrix[,match(sample,colnames(matrix))]
matrix_edge <- tissue_matrix[rowSums(tissue_matrix)>0,]
ref <- reference[match(colnames(matrix_edge),reference[,1]),]
Meta_data_sort <- ref[order(ref$Predict_Tissue),]


fdrThres=0.05

####exclude tissues with the same biological background(such as when you do on small intestine, you have to remove all other gut tissues)
select_table <- group[group[,1] == tissue,]
tissue_table <- group[group[,1] == tissue |group[,2] == tissue,]
exclude_tissue <- unique(
  group[group[, 2] %in% select_table[, 2], 2]
)
exclude_tissue <- setdiff(exclude_tissue, select_table[, 2])
temp_meta <- ref[!ref$Predict_Tissue %in% exclude_tissue, ]
temp_meta$conditions <- ifelse(
  temp_meta$Predict_Tissue %in% select_table[, 2],
  "Tissue",
  "Others"
)
conditions <- factor(temp_meta$conditions)
temp_matrix <- matrix_edge[, temp_meta$BioSample, drop = FALSE]


####wilcoxon DEG  ##https://link.springer.com/article/10.1186/s13059-022-02648-4
y <- DGEList(counts=temp_matrix,group=conditions)
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)
pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
outRst <- outRst[order(outRst$FDR,decreasing = F),]
output_dir <- "./tissue_specific_gene_all/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.table(outRst, file=paste0("./tissue_specific_gene_all/",tissue_name,"_Wilcoxon.rst.tsv"),sep="\t", quote=F,row.names = T,col.names = T)
