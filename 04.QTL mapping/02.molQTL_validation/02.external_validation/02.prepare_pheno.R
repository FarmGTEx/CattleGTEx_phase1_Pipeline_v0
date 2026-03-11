library(data.table)
library(edgeR)
library(data.table)


"%&%" = function(a, b) { paste0(a, b) }
inverse_normal_transform = function(x) {
    qnorm(rank(x) / (length(x)+1))
}


TPM <- readRDS("/faststorage/project/farmgtex/pipeline/QTL_validation/TPM_matrix.rds")
count <- readRDS("/faststorage/project/farmgtex/pipeline/QTL_validation/count_matrix.rds")
leave <- fread("/faststorage/project/farmgtex/pipeline/QTL_validation2/overlap_samples.txt",header = F,data.table = F)[,1]
tss_annot_file <- "/faststorage/project/farmgtex/pipeline/eQTL/TSS.gtf"
gene_filter <- "/faststorage/project/farmgtex/pipeline/genechosen"

colnames(TPM) <- gsub("_batch1$", "", colnames(TPM))
TPM <- TPM[,colnames(TPM) %in% leave]

colnames(count) <- gsub("_batch1$", "", colnames(count))
Counts <- count[,colnames(count) %in% leave]

filter = fread(gene_filter,header = F,data.table = F)
samids = colnames(Counts) # sample id
expr_counts = Counts[rownames(Counts) %in% filter[,1],]
expr = DGEList(counts=expr_counts) # counts
nsamples = length(samids) # sample number
ngenes = nrow(expr_counts) # gene number



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
TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = inverse_normal_transform)) 

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
fwrite(bed, "all.expr_tmm_inv.bed", sep = "\t")
system("bgzip -f all.expr_tmm_inv.bed")
system("tabix -p bed all.expr_tmm_inv.bed.gz")