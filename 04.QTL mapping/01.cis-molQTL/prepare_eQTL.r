options(stringsAsFactors = FALSE)
ARGS <- commandArgs(trailingOnly = TRUE)
library(edgeR)
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
file_tpm <- "/faststorage/project/farmgtex/pipeline/eQTL/" %&% scenario_select %&% "/TPM_QTL.rds"
file_counts <- "/faststorage/project/farmgtex/pipeline/eQTL/" %&% scenario_select %&% "/count_QTL.rds"
tss_annot_file <- "/faststorage/project/farmgtex/pipeline/eQTL/TSS.gtf"
gene_filter <- "/faststorage/project/farmgtex/pipeline/genechosen"
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
### main program
#----------------------------------------------------------------------------
# Input data

##TPM matrix. Row is transcript, column is sample; rowname is gene id, colname is sample id
TPM = readRDS(file_tpm)
## Read counts matrix. Row is gene, column is sample; rowname is gene id, colname is sample id
Counts = readRDS(file_counts)
## Gene filter. Only consider lncRNA and PCG in this step
filter = fread(gene_filter,header = F,data.table = F)

################10.2
## 1. prepare TMM
samids = colnames(Counts) # sample id
expr_counts = Counts[rownames(Counts) %in% filter[,1],]
expr = DGEList(counts=expr_counts) # counts
nsamples = length(samids) # sample number
ngenes = nrow(expr_counts) # gene number


# calculate TMM
y = calcNormFactors(expr, method="TMM")
TMM = cpm(y,normalized.lib.sizes=T)

# expression thresholds

count_threshold = 6
tpm_threshold = 0.1
sample_frac_threshold = 0.2
sample_count_threshold = 10


#keep the gene with >=0.1 tpm and >=6 read counts in >=20% samples.
expr_tpm = TPM[rownames(expr_counts),samids]
tpm_th = rowSums(expr_tpm >= tpm_threshold)
count_th = rowSums(expr_counts >= count_threshold)
ctrl1 = tpm_th >= (sample_frac_threshold * nsamples)
ctrl2 = count_th >= (sample_frac_threshold * nsamples)
mask = ctrl1 & ctrl2
TMM_pass = TMM[mask,] ##row is gene; column is sample

###expression values (TMM) were inverse normal transformed across samples.
TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = inverse_normal_transform)) #apply to each row, each row represents one gene, observed values for all the samples. scale across samples.
#----------------------------------------------------------------------------
### 2. prepare bed file
dir.create("bed",showWarnings=F)

region_annot = fread(tss_annot_file) # load gtf file
geneid = region_annot$gene_id

expr_matrix = TMM_inv[rownames(TMM_inv) %in% geneid,] # expr_matrix TMM_inv

# prepare bed file for tensorQTL
bed_annot = region_annot[region_annot$gene_id %in% rownames(expr_matrix),]
bed = data.frame(bed_annot,expr_matrix[bed_annot$gene_id,],check.names = FALSE)
bed = bed[bed[,1] %in% as.character(1:29),]
bed[,1] = as.numeric(bed[,1])
bed = bed[order(bed[,1],bed[,2]),]
colnames(bed)[1] = "#Chr"

# output bed file
fwrite(bed,file = "./bed/rna.expr_tmm_inv.bed", sep = "\t")
system("bgzip -f ./bed/rna.expr_tmm_inv.bed")
system("tabix -p bed ./bed/rna.expr_tmm_inv.bed.gz")
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
cat("done.\n")
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
