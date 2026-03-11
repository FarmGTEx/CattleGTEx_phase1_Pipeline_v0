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
file_counts <- "/faststorage/project/farmgtex/pipeline/enQTL/" %&% scenario_select %&% "/count_QTL.rds"
tss_annot_file <- "/faststorage/project/farmgtex/pipeline/enQTL/Enhancer_pos.gtf"
gene_filter <- "/faststorage/project/farmgtex/pipeline/enQTL/enhancer.bed"
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
### main program
#----------------------------------------------------------------------------
# Input data


## Read counts matrix. Row is gene, column is sample; rowname is gene id, colname is sample id
Counts = readRDS(file_counts)


################10.2
## 1. prepare TMM
samids = colnames(Counts) # sample id
expr_counts = Counts
expr = DGEList(counts=expr_counts) # counts
nsamples = length(samids) # sample number
ngenes = nrow(expr_counts) # gene number


# calculate TMM
y = calcNormFactors(expr, method="TMM")
TMM = cpm(y,normalized.lib.sizes=T)

# expression thresholds

count_threshold = 6
sample_frac_threshold = 0.2
sample_count_threshold = 10


#keep the gene with >=0.1 tpm and >=6 read counts in >=20% samples.
count_th = rowSums(expr_counts >= count_threshold)
ctrl = count_th >= (sample_frac_threshold * nsamples)
mask = ctrl
TMM_pass = TMM[mask,] ##row is gene; column is sample

###expression values (TMM) were inverse normal transformed across samples.
TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = inverse_normal_transform)) #apply to each row, each row represents one gene, observed values for all the samples. scale across samples.
#----------------------------------------------------------------------------
### 2. prepare bed file
dir.create("bed",showWarnings=F)

region_annot = fread(tss_annot_file,data.table = F) # load gtf file
gene_filter = fread(gene_filter,data.table = F)
region_annot = region_annot[region_annot[,4] %in% gene_filter[,4],]
enhancer_id = region_annot$enhancer_id
expr_matrix = TMM_inv[rownames(TMM_inv) %in% enhancer_id,] # expr_matrix TMM_inv
bed_annot <- region_annot[region_annot$enhancer_id %in%rownames(expr_matrix),]
# prepare bed file for tensorQTL
bed = data.frame(bed_annot,expr_matrix[bed_annot$enhancer_id,],check.names = FALSE)
bed = bed[bed[,1] %in% as.character(1:29),]
bed[,1] = as.numeric(bed[,1])
bed = bed[order(bed[,1],bed[,2]),]
colnames(bed)[1] = "#Chr"
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
