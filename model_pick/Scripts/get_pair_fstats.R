library(sgeostat)
library(pbkrtest)
library(mvoutlier)
library(boot)

args = commandArgs(trailingOnly=TRUE)
print(args)
p_edges = args[1]
p_all_models = args[2]
assign("pval_thr", as.numeric(args[3]), envir = .GlobalEnv)
out = args[4]

# p_edges = 'edges_0.2_0.txt'
# p_all_models = 'all_models.txt'
# pval_thr = 0.001
# out = './'

pull_edges <- function(p_edges, p_all_models){
  edges = read.table(p_edges, stringsAsFactors = F)
  all_models = read.table(p_all_models, header = 1, stringsAsFactors = F)
  pair_models = data.frame(matrix(ncol = 6, nrow = 0))
  colnames(pair_models) = c('Predictor', 'Response', 'F_stat', 'Coef_p', 'Intercept', 'P_value')
  for (i in rownames(edges)){
    pred_resp = list(
      direct = c(edges[i,1],edges[i,2]),
      reverse = c(edges[i,2],edges[i,1])
    )
    for (k in c('direct', 'reverse')) {
      predictor = pred_resp[[k]][1]
      response = pred_resp[[k]][2]
      temp = all_models[all_models['Predictor'] == predictor,]
      if ( nrow(temp) == 0 ) {
        next
      }
      temp = temp[temp['Response'] == response,]
      if ( nrow(temp) == 0) {
        next
      }
      if (temp$P_value < pval_thr) {
        pair_models = rbind(pair_models,temp[1,], make.row.names = F)
      }
    }
  }
  if (nrow(pair_models) == 0){
    return('fail')
  }
  return(pair_models)
}

# Case of empty pair correlation graph
if (readLines(p_edges, n=1) == "No edges based on pair correlations"){
  outfile = paste0(out,'/pair_models_', toString(pval_thr), '.txt')
  writeLines('No significant models', outfile)
} else {
# Case of non-empty pair correlation graph
  
  fstats = pull_edges(p_edges, p_all_models)
  
  # This line used to work, but now it doesn't.
  # hmmmmmmmmm
  # mystery
  # outfile = paste0(out,'/pair_models_', toString(fstat_thr), '.txt')
  outfile = paste0('pair_models_', toString(pval_thr), '.txt')
  setwd(out)
  print(file.exists(out))
  if (fstats == 'fail') {
    writeLines('No significant models', outfile)
  } else {
    write.table(fstats, outfile, sep='\t', row.names = F, quote = F)
  }
}
