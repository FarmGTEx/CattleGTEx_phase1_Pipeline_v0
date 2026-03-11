library(data.table)
library(dplyr)
"%&%" = function(a, b) { paste0(a, b) }
ARGS <- commandArgs(trailingOnly = TRUE)
scenario_select <- ARGS[1] 


###using Leafcutter to generate splicing ratio matrix
setwd(paste0("/faststorage/project/farmgtex/pipeline/sQTL/",scenario_select))
dir.create("bed",showWarnings=F)
clusters <- fread("cluster_list.txt",header = F, data.table = F)[,1]
samples <- fread("sample_list.txt",header = F, data.table = F)[,1]
all_samples <- fread("/faststorage/project/farmgtex/pipeline/sQTL/IDs",header = F, data.table = F)[,1]
#system("awk 'NR==FNR {clusters[$1]; next} $1 in clusters' cluster_list.txt /faststorage/project/farmgtex/pipeline/splicing/cattle_30_perind.counts > test")
sub <- read.table("test", header = FALSE)
colnames(sub)[2:ncol(sub)] <- all_samples
colnames(sub)[1] <- "chrom"
sample_select <- cbind(sub$chrom,sub[,colnames(sub)%in% samples])
colnames(sample_select)[1] <- "chrom"
write.table(sample_select,file = "filtered_matrix.bed",col.names= T,row.names = F,quote = F, sep = " ")
system("gzip -f filtered_matrix.bed")
system("python /faststorage/project/farmgtex/gtex/tools/leafcutter/scripts/prepare_phenotype_table.py ./filtered_matrix.bed.gz")
system("rm temp.matrix")
system("ls -v filtered_matrix.bed.gz.qqnorm_* | grep -E '_[0-9]+$' | while read chr; do cat $chr >> temp.matrix; done")
system("sed '0,/^#/!{/^#/d}' temp.matrix | grep -v \"Chr\" | sort -V -k1,1 -k2,2n -k3,3n > sorted.bed")
bed <- fread("sorted.bed",header = F, data.table = F)
annotation <- fread("/faststorage/project/farmgtex/pipeline/eQTL/TSS.gtf",header = T,data.table =F)
group_file <- fread("cluster_groups.txt",header = F,data.table =F)
info <- bed[,1:4]
info$gene_id <- group_file[match(info[,4],group_file[,1]),2]
bed[,2] <- annotation[match(info$gene_id,annotation$gene_id),2]
bed[,3] <- annotation[match(info$gene_id,annotation$gene_id),3]
bed_clean <- na.omit(bed)
#group_file <- group_file[match(bed_clean[,4,],group_file[,1]),]
group_file <- group_file[group_file[,1] %in% bed_clean[,4,],]
bed_sorted <- bed_clean[order(bed_clean$V1, bed_clean$V2, bed_clean$V3), ]
group_file_sorted<-group_file[match(bed_sorted[,4],group_file[,1]),]
group_list <- unique(group_file_sorted[, 2])
bed_sorted$Group <- group_file_sorted[, 2]
bed_sorted$Group <- factor(bed_sorted$Group, levels = group_list)
bed_sorted_final <- bed_sorted[order(bed_sorted$Group), ]
group_info <- bed_sorted_final[,c(4,ncol(bed_sorted_final))]
bed_sorted_final$Group <- NULL
write.table(bed_sorted_final,"./bed/"%&% scenario_select %&%"_sorted.bed",quote = F,col.names = F, row.names = F, sep= "\t")
write.table(group_info,"./sample_group.txt",quote = F,col.names = F, row.names = F, sep= "\t")
system("less temp.matrix | grep \"Chr\" | uniq  > ./bed/header")
system("cat ./bed/header ./bed/"%&% scenario_select %&%"_sorted.bed > ./bed/"%&% scenario_select %&%"_filted_qqnorm.bed")
system("rm ./bed/header ./bed/"%&% scenario_select %&%"_sorted.bed" )
system("bgzip -f ./bed/"%&% scenario_select %&%"_filted_qqnorm.bed")
system("tabix -f -p bed ./bed/"%&% scenario_select %&%"_filted_qqnorm.bed.gz")