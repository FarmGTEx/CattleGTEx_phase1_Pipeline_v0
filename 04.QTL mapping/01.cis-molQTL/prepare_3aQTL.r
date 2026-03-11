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
file_APA <- "/faststorage/project/farmgtex/pipeline/polyA/" %&% scenario_select %&% "/" %&%scenario_select %&%"_final_predict_clean.txt"
tss_annot_file <- "/faststorage/project/farmgtex/pipeline/3aQTL/" %&% scenario_select %&% "/annotation.gtf"
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
### main program
#----------------------------------------------------------------------------
# Input data

## Read counts matrix. Row is gene, column is sample; rowname is gene id, colname is sample id
APA = fread(file_APA,header = T, data.table = F)
rownames(APA) <- APA[,1]
APA <- APA[,-1]


################10.2
## 1. prepare TMM
samids = colnames(APA) # sample id
nsamples = length(samids) # sample number
ngenes = nrow(APA) # gene number


###expression values (TMM) were inverse normal transformed across samples.
APA_inv = t(apply(APA, MARGIN = 1, FUN = inverse_normal_transform)) #apply to each row, each row represents one gene, observed values for all the samples. scale across samples.
rows_to_keep <- apply(APA_inv, 1, function(x) length(unique(x)) > 1)
APA_inv_filtered <- APA_inv[rows_to_keep, ]
#----------------------------------------------------------------------------
### 2. prepare bed file
dir.create("bed",showWarnings=F)

bed_annot = fread(tss_annot_file) # load gtf file
featureid = bed_annot$PolyAloci
bed_annot <- bed_annot[bed_annot$PolyAloci %in% rownames(APA_inv_filtered),]
bed = data.frame(bed_annot,APA_inv_filtered[bed_annot$PolyAloci,],check.names = FALSE)
bed = bed[bed[,1] %in% as.character(1:29),]
bed[,1] = as.numeric(bed[,1])
bed = bed[order(bed[,1],bed[,2]),]
colnames(bed)[1] = "#Chr"
bed <- bed[apply(bed[, 5:ncol(bed)], 1, function(x) {
  any(x != 0)
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
