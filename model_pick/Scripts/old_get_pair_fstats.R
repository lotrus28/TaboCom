library(sgeostat)
library(pbkrtest)
library(mvoutlier)
library(boot)

args = commandArgs(trailingOnly=TRUE)
print(args)

edges = args[1]
counts = args[2]
tax_code = args[3]
out = args[4]

assign("fstat_thr", as.numeric(args[5]), envir = .GlobalEnv)

edges_to_sig_cor = function(p_edges) {
  edges = read.csv(p_edges, sep = '\t',stringsAsFactors = F)
  sig.cor = data.frame()
  for (i in 1:nrow(edges)) {
    sig.cor[edges[i,1],edges[i,2]] = edges[i,3]
    sig.cor[edges[i,2],edges[i,1]] = edges[i,3]
  }
  sig.cor[is.na(sig.cor)] <- 0
  return(sig.cor)
}

get_pair_models_table = function(path_to_cor = NA, path_to_counts = NA, 
                                  cor_table = NA, counts = NA, bootnum, boot_pV) {
  
  # in: sig.cor table from "get_sig_cor"
  # out: list of sign-ly correlated organisms
  # (every odd element is #1 in a pair, every even -- #2)
  get_sig_pairs <- function(sig_cor){
    # Shouldn't sig_cor be a square matrix?
    L = ncol(sig_cor)
    r_c_names = rownames(sig_cor)
    sig_pairs = c()
    # Go othrough matrix of sig_corrs diagonally
    # and remember significant pairs
    
    # In case something goes wrong: use which instead
    # sig_cor[lower.tri(sig_cor)] = 0
    # indices = which(sig_cor != 0, arr.ind=TRUE)
    # rownames(indices) = c(1:nrow(indices))
    
    for (i in 1:(L-1))
    {
      for (j in (1+i):L)
      {
        if (sig_cor[r_c_names[i],r_c_names[j]] != 0)
        {
          new_pair = c(r_c_names[i], r_c_names[j])
          sig_pairs = c(sig_pairs, new_pair)
        }
      }
    }
    return(sig_pairs)
  }
  
  # in: list of pairs from sig_pairs, count data,
  # pV to assess model quality, # of bootstrap iterations (pV >= 1/num_boot)
  # out: a table of non-categorial pair-models
  # (Each cell has a list w. 1 element ? lm object)
  # (Row/Colnames are taken from sig_pairs - rows = predictors, cols=response)
  produce_pair_models <- function(sig_pairs, initial_data, pV_cutoff, num_boot) {
    
    get_bootstrap <- function(d,i){
      
      sim.data <- function(data,mle){
        
        rden <- function(n, den) {
          # diffs <- diff(den$x)
          # Making sure we have equal increments
          # stopifnot(all(abs(diff(den$x) - mean(diff(den$x))) < 1e-9))
          total <- sum(den$y)
          den$y <- den$y / total
          ydistr <- cumsum(den$y)
          yunif <- runif(n)
          indices <- sapply(yunif, function(y) min(which(ydistr > y)))
          x <- den$x[indices]
          
          return(x)
        }
        resp = as.vector(data[,2])
        pred = as.vector(data[,1])
        den_resp = density(resp, from = 0)
        den_pred = density(pred, from = 0)
        N = nrow(data)
        sim_resp = rden(N,den_resp)
        sim_pred = rden(N,den_pred)
        temp = as.data.frame(cbind(sim_resp, sim_pred))
        return(temp)
      }
      
      stat<- function(data) {
        model <- lm(data[,1] ~ data[,2], data)
        F_stat <- summary(model)$fstatistic
        if(is.null(F_stat)) {
          F_score <- 1
        } else {
          F_score <- pf(F_stat[1],F_stat[2],F_stat[3],lower.tail = FALSE)
        }
        return(F_score)
      }
      
      temp <- boot(d, statistic=stat, R=i,
                   ran.gen=sim.data,
                   sim="parametric", mle = 1)
      F_scores <- temp$t
      
      return(F_scores)
      
    }

    pred_names <- unique(sig_pairs[seq(1, length(sig_pairs), 2)])
    res_names <- unique(sig_pairs[seq(2, length(sig_pairs), 2)])
    
    # Empty matrix to contain linear models
    all_models = matrix(list(NA), ncol = length(res_names), nrow = length(pred_names))
    
    rownames(all_models) <- pred_names
    colnames(all_models) <- res_names
    L = length(sig_pairs)/2.
    count_models = 0
    
    # Fill the matrix with models
      for (i in 1:L) 
    {
      # Reverse models are not checked
      predictor = sig_pairs[(i*2)-1]
      response = sig_pairs[i*2]
      
      temp = as.data.frame(t(rbind(initial_data[response,],initial_data[predictor,])))
      colnames(temp) <- c('resp','pred')
      
      # Remove all zero-pairs but one
      temp_not_zeros = temp[apply(temp,1,function(row) !all(row ==0 )),]
      temp_zeros = temp[apply(temp,1,function(row) all(row ==0 )),]
      if (nrow(temp_zeros) > 1) {
        temp = rbind(temp_not_zeros, c(0,0))
      } else {
        temp = rbind(temp_not_zeros, temp_zeros)
      }
      
      # Get rid of outliers
      # If insufficient data - go to next pair
      outs <- tryCatch({
        dd.plot(temp)$outliers
      }, error = function(e) {
		print(e)
        -1
      })

      if (outs == -1) next
      temp <- temp[!outs,]
      
      # Skip if a column is full of zeroes
      if (length(temp[,colSums(temp^2) == 0]) > 0) next
      text = paste('resp', 'pred', sep = ' ~ ')
      param <- list(formula = text, data = temp)
      model <- do.call(lm, param)
      # f - how much a model with regression is better than
      # an intercept only model
      tryCatch({
        f = summary(model)$fstatistic
      }, warning=function(w) {
        message("Warning: ", conditionMessage(w))
        print(text)
      })
      
      # Calculate pValue of a model
      # based on the # of bootstrapped models
      # w. better F-statistics
      sim_F <- get_bootstrap(temp,num_boot)
      true_F = pf(f[1],f[2],f[3],lower.tail = FALSE)
      pV = length(sim_F[true_F > sim_F]) / length(sim_F)
      # print(paste0('pValue: ', pV))
      
      if (pV< fstat_thr) 
      {
        count_models = count_models+1
        # Somehow predict wont work if I touch $call
        # model$call <- NULL
        
        # At first I used this way of saving lists to a dataframe
        # all_models[[predictor,response]] <-  list(c(model))
        # But it is unstable. Why so? No idea
        
        # So I decided to make a list-containing matrix
        # and fill it with lists
        all_models[predictor,response] <-  list(model)
      }
    }
    
    print(paste0('Models built: ', count_models))
    # Remove empty columns and rows
    all_models = all_models[ , ! apply( all_models , 2 , function(x) all(is.na(x)) ) ]
    all_models = all_models[! apply( all_models , 1 , function(x) all(is.na(x)) ) ,]
    return(as.data.frame(all_models))
  }
  
  if (!(is.na(path_to_cor))) {
    sig_cor <- read.table(path_to_cor,sep="\t", header=TRUE, row.names = 1)
  } else {
    sig_cor <- cor_table
  }
  if (!(is.na(path_to_counts))) {
    counts <-  read.table(path_to_counts, sep = '\t', header=TRUE, row.names = 1)
  } else {
    counts <- counts
  }
  
  sig_pairs  <-  get_sig_pairs(sig_cor)
  # sig_pairs = get_sig_sets(sig_cor,2)
  # In pairs for loop calculation is faster
  all_models <-  produce_pair_models(sig_pairs, counts, boot_pV, bootnum)
  return(all_models)
}

get_pair_fstats = function(models) {
  
  meta_models = data.frame(Predictor = 'pred', Response = 'resp', F_Stat = 1, Coef_p = 0, Intercept = 0, stringsAsFactors = FALSE)
  
  for (i in 1:dim(models)[1]) {
    
    for (j in 1:dim(models)[2]){
      
	  if (is.null(models[i,j][[1]])) {
          models[i,j][[1]] <- NA
        }
	  
      temp = NA
      if (!(is.na(models[i,j][[1]]))) {
        f = summary(models[i,j][[1]])$fstatistic
        pF = pf(f[1],f[2],f[3],lower.tail = FALSE)
        coef_p = coef(models[i,j][[1]])[2]
        intercept = coef(models[i,j][[1]])[1]
        temp = c(dimnames(models)[[1]][i],dimnames(models)[[2]][j],pF, coef_p, intercept)
      } 
      if (!is.na(temp)) meta_models <- rbind(meta_models, temp)
    }
  }
  meta_models <- meta_models[-1,]
  rownames(meta_models) <- 1:nrow(meta_models) 
  
  return(meta_models)
}


# Case of empty pair correlation graph
if (readLines(edges, n=1) == "No edges based on pair correlations"){
  outfile = paste0(out,'/pair_models_', toString(fstat_thr), '.txt')
  writeLines('No significant models', outfile)
} else {
# Case of non-empty pair correlation graph
  sig.cor = edges_to_sig_cor(edges)
  
  # Taxonomies to IDs
  spc = read.table(tax_code , stringsAsFactors = F, row.names = 2)
  colnames(sig.cor) = spc[colnames(sig.cor),]
  rownames(sig.cor) = spc[rownames(sig.cor),]
  sig.cor = sig.cor[order(rownames(sig.cor)) , order(colnames(sig.cor))]
  
  counts = read.table(counts, stringsAsFactors = F, row.names = 1)
  
  
  fstats = tryCatch({
    pair_mod = get_pair_models_table(cor_table = sig.cor, counts = counts, bootnum = 1/fstat_thr, boot_pV = fstat_thr)
    # IDs back to taxonomies
    spc = read.table(tax_code, stringsAsFactors = F, row.names = 1)
    colnames(pair_mod) = spc[colnames(pair_mod),]
    rownames(pair_mod) = spc[rownames(pair_mod),]
    fstats = get_pair_fstats(pair_mod)
  }, error = function(err) {
    return('fail')
  })
  
  # This line used to work, but now it doesn't.
  # hmmmmmmmmm
  # mystery
  # outfile = paste0(out,'/pair_models_', toString(fstat_thr), '.txt')
  outfile = paste0('pair_models_', toString(fstat_thr), '.txt')
  setwd(args[4])
  print(file.exists(args[4]))
  if (fstats == 'fail') {
    writeLines('No significant models', outfile)
  } else {
    write.table(fstats, outfile, sep='\t', row.names = F, quote = F)
  }
}
