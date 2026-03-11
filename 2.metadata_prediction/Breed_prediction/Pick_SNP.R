args <- commandArgs(trailingOnly = TRUE)
N=args[1]
library(data.table)
info <- fread("FST_AED_PLSR.txt", header = F, data.table = F)
colnames(info) <- c("SNP","FST","AED","Pos","PLSR")
AED <- info[order(info$AED,decreasing =T),]
AED <- AED[1:N,]
FST <- info[order(info$FST,decreasing =T),]
FST <- FST[1:N,]
PLSR <- info[order(info$PLSR,decreasing =T),]
PLSR <- PLSR[1:N,]
write.table(AED$SNP,paste0("AED_",N,"snp.txt"),col.names = F, row.names= F , quote = F , sep = "\t")
write.table(FST$SNP,paste0("FST_",N,"snp.txt"),col.names = F, row.names= F , quote = F , sep = "\t")
write.table(PLSR$SNP,paste0("PLSR_",N,"snp.txt"),col.names = F, row.names= F , quote = F , sep = "\t")


system(paste0("plink --bfile impute_samples --recode A --extract FST_",  N, "snp.txt --threads 12 --cow --out FST_",  N, "_ped"))
system(paste0("plink --bfile impute_samples --recode A --extract AED_",  N, "snp.txt --threads 12 --cow --out AED_",  N, "_ped"))
system(paste0("plink --bfile impute_samples --recode A --extract PLSR_", N, "snp.txt --threads 12 --cow --out PLSR_", N, "_ped"))



