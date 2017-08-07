args = commandArgs(trailingOnly=TRUE)

print(args)

# args = c('pvals_two_sided.txt', 'cor_mat_SparCC.out', 'tax_code.txt', "E://Lab//FHM//AGP//Parallel//py", '0.01')

assign("pval_thr", as.numeric(args[5]), envir = .GlobalEnv)

get_sig_cor = function(p_pval = NA, p_corr = NA, val = NA, corr = NA, thr_pval = pval_thr){
  if (!(is.na(p_pval))) {
    pvals=read.table(p_pval,header=TRUE, row.names = 1, sep="\t")
    cor = read.table(p_corr,header=TRUE, row.names = 1, sep="\t")
  }
  else {
    pvals = val
    cor = corr
  }
  pvals.mat=pvals[,1:ncol(pvals)]
  pvals.mat[pvals.mat==0]=0.000000001
  # convert into significance
  sig.mat= -1*log10(pvals.mat)
  
  # thr_pval = (-1*log10(2./(ncol(pvals)**2)))/7
  # thr_pval = -1*log10(thr_pval)
  
  # remove all edges with significance below threshold
  adj.mat=sig.mat
  adj.mat[adj.mat<thr_pval]=0
  adj.mat[adj.mat>=thr_pval]=1
  adj.mat=as.matrix(adj.mat)
  # adj.mat=adj.mat[,colSums(adj.mat^2) !=0]
  
  sign.cor=as.matrix(cor)
  for (i in 1:nrow(adj.mat)) {
    for (j in 1:ncol(adj.mat)) {
      sign.cor[i,j] <- adj.mat[i,j]*sign.cor[i,j]
    }
  }
  return(sign.cor)
}


if (length(args)==5) {
	sign.cor = get_sig_cor(p_pval =  args[1], p_corr =  args[2], thr_pval = as.numeric(args[5]))
}

spc = read.table(args[3], stringsAsFactors = F, row.names = 1)
rownames(sign.cor)  <- spc[rownames(sign.cor),]
colnames(sign.cor)  <- spc[colnames(sign.cor),]

# write.table(sc, paste0(args[4],'/sig_cor.txt'), sep='\t', col.names=NA, quote = F)
outfile = paste0('sig_cor_', toString(pval_thr), '.txt')
setwd(args[4])
write.table(sign.cor, outfile, sep='\t', col.names=NA, quote = F)
