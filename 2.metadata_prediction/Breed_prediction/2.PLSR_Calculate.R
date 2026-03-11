library(pls)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

file=args[1]
library(data.table)
matrix <- fread(paste0("/temp/",file,".raw"), header =T , data.table = F)#convert the file to PLINK .raw format (additive genotype coding) 
matrix = matrix[,-c(2:6)]
rownames(matrix) <- matrix[,1]
matrix = matrix[,-1]

#match breed info to samples 
RNA_sample <- fread("")  #metadata of RNA-seq samples
WGS_sample <- fread("")  #metadata of WGS samples from reference panel


colnames(RNA_sample) <- colnames(WGS_sample)
all_sample <- rbind(RNA_sample,WGS_sample)
matrix$breed <- all_sample[match(rownames(matrix),all_sample[,1]),2]
matrix$breed<-as.factor(matrix$breed)
matrix$breed <- as.numeric(matrix$breed)
matrix <- as.data.frame(sapply(matrix, as.numeric))

#fit a PLSR model to estimate coefficients for each SNP and save the results
plsr_model <- plsr(breed ~ ., data = matrix, ncomp = 100)
loadings <- loadings(plsr_model)
ncomp <- ncol(loadings)
VIP <- numeric(length = nrow(loadings))
for (j in 1:nrow(loadings)) {
  VIP[j] <- sum((loadings[j, ]^2))
}

write.table(cbind(colnames(matrix[,1:ncol(matrix)-1]),VIP),paste0(file,".txt"),quote = F,sep = "\t",row.names = F,col.names = F)