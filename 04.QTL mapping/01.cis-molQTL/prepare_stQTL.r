options(stringsAsFactors = FALSE)
ARGS <- commandArgs(trailingOnly = TRUE)
library(data.table)

#----------------------------------------------------------------------------
### functions
"%&%" = function(a, b) { paste0(a, b) }
inverse_normal_transform = function(x) {
    qnorm(rank(x) / (length(x)+1))
}
#----------------------------------------------------------------------------
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
### data input <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ARGS <- commandArgs(trailingOnly = TRUE)

scenario_select <- ARGS[1] # Prefix of output file, like tissue name  
file_ratio <- "rna_ratio_predict.txt"
tss_annot_file <- "/faststorage/project/farmgtex/pipeline/eQTL/TSS.gtf"
gene_filter <- "/faststorage/project/farmgtex/pipeline/genechosen"
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
### main program
#----------------------------------------------------------------------------
# Input data
## Read counts matrix.

filter = fread(gene_filter,header = F,data.table = F)

ratio_matrix <- fread(file_ratio,header = T,data.table = F)
rownames(ratio_matrix) <- ratio_matrix[,1]
ratio_matrix <- ratio_matrix[,-1]
ratio_matrix <- ratio_matrix[rownames(ratio_matrix) %in% filter[,1],]

###expression values (TMM) were inverse normal transformed across samples.
ratio_inv = t(apply(ratio_matrix, MARGIN = 1, FUN = inverse_normal_transform)) #apply to each row, each row represents one gene, observed values for all the samples. scale across samples.
#----------------------------------------------------------------------------
### 2. prepare bed file
dir.create("bed",showWarnings=F)

region_annot = fread(tss_annot_file) # load gtf file
gene_id = region_annot$gene_id

ratio_matrix = ratio_inv[rownames(ratio_inv) %in% gene_id,] # expr_matrix TMM_inv
ratio_matrix <- ratio_matrix[rowSums(ratio_matrix != 0) > 0, ]
# prepare bed file for tensorQTL
bed_annot = region_annot[region_annot$gene_id %in% rownames(ratio_matrix),]
bed = data.frame(bed_annot,ratio_matrix[bed_annot$gene_id,],check.names = FALSE)
bed = bed[bed[,1] %in% as.character(1:29),]
bed[,1] = as.numeric(bed[,1])
bed = bed[order(bed[,1],bed[,2]),]
colnames(bed)[1] = "#Chr"
bed <- bed[apply(bed[, 5:ncol(bed)], 1, function(x) {
  freq <- table(x)
  !(length(freq) == 2 && any(freq == 1))
}), ]
# output bed file
fwrite(bed,file = "./bed/" %&% scenario_select %&% ".expr_tmm_inv.bed", sep = "\t")
system("bgzip -f ./bed/" %&% scenario_select %&% ".expr_tmm_inv.bed")
system("tabix -p bed ./bed/" %&% scenario_select %&% ".expr_tmm_inv.bed.gz")
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
cat("done.\n")
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
