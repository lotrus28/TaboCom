library(DESeq)

args = commandArgs(trailingOnly=TRUE)
assign("path_to_counts", args[1], envir = .GlobalEnv)
assign("wd", args[2], envir = .GlobalEnv)

setwd(wd)

data = read.table(path_to_counts, header=T, sep = '\t', row.names = 1)
data[is.na(data)] <- 0
conds = rep('cond1',ncol(data))

otu = newCountDataSet(data, conds)
otu <- estimateSizeFactors( otu )
otu <- estimateDispersions( otu, fitType="local")

normalizedCounts <- t( t(counts(otu)) / sizeFactors(otu) )

write.table(normalizedCounts, file = 'train_norm.txt', sep = '\t', col.names = NA, quote=F)
