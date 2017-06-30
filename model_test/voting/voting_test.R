# args = commandArgs(trailingOnly=TRUE)
args = c('CROSS_PARAMETERS.txt',
         './FINAL_OTU_TABLE_ibd/',
         './FINAL_OTU_TABLE_heal/',
         './',
         'test')

params = args[1]
I_folders = args[2]
H_folders = args[3]
out = args[4]
test_or_valid = args[5]

parse_info <- function(row, I_folders, H_folders, test_or_valid){
  
  common_tax_code <- function(I_tax_code,H_tax_code){
    
    c_tax = data.frame(NA, NA , intersect(I_tax_code[[2]],H_tax_code[[2]]), NA)
    colnames(c_tax) = c('H_ID','I_ID', 'Taxon', 'Common')
    
    # Add common taxons
    for (i in rownames(c_tax)){
      taxon = c_tax[i,'Taxon']
      c_tax[i,'I_ID'] = I_tax_code[which(I_tax_code[[2]] == taxon),1]
      c_tax[i,'H_ID'] = H_tax_code[which(H_tax_code[[2]] == taxon),1]
    }
    
    # Add H-unique taxons
    h_tax = data.frame(NA, NA , setdiff(H_tax_code[[2]],I_tax_code[[2]]), NA)
    colnames(h_tax) = c('H_ID','I_ID', 'Taxon', 'Common')
    for (i in rownames(h_tax)){
      taxon = h_tax[i,'Taxon']
      h_tax[i,'H_ID'] = H_tax_code[which(H_tax_code[[2]] == taxon),1]
    }
    
    # Add I-unique taxons
    i_tax = data.frame(NA, NA , setdiff(I_tax_code[[2]],H_tax_code[[2]]), NA)
    colnames(i_tax) = c('H_ID','I_ID', 'Taxon', 'Common')
    for (i in rownames(i_tax)){
      taxon = i_tax[i,'Taxon']
      i_tax[i,'I_ID'] = I_tax_code[which(I_tax_code[[2]] == taxon),1]
    }
    
    # Merge all tables
    c_tax = rbind(c_tax, h_tax, i_tax)
    
    # Add common tax_ID
    c_tax['Common'] = paste0('tax', rownames(c_tax))
    
    return(c_tax)
  }
  
  # Find corresponding models
  if (length(colnames(row)) == 13) {
    H_Taxon_trim = row$H_Taxon_trim
    H_SparCC_pval = row$H_SparCC_pval
    H_SparCC_cor = row$H_SparCC_cor
    H_Tax_adjacency = row$H_Tax_adjacency
    H_Pair_fstat = row$H_Pair_fstat 
    I_Taxon_trim = row$I_Taxon_trim
    I_SparCC_pval = row$I_SparCC_pval
    I_SparCC_cor = row$I_SparCC_cor
    I_Tax_adjacency = row$I_Tax_adjacency
    I_Pair_fstat = row$I_Pair_fstat
    RMSE_sign = row$RMSE_sign
  } else {
    H_Taxon_trim = row$Taxon_trim
    H_SparCC_pval = row$SparCC_pval
    H_SparCC_cor = row$SparCC_cor
    H_Tax_adjacency = row$Tax_adjacency
    H_Pair_fstat = row$Pair_fstat 
    I_Taxon_trim = row$Taxon_trim
    I_SparCC_pval = row$SparCC_pval
    I_SparCC_cor = row$I_SparCC_cor
    I_Tax_adjacency = row$Tax_adjacency
    I_Pair_fstat = row$Pair_fstat
    RMSE_sign = row$RMSE_sign
  }
  
  for (condition in c('H','I')){
    
    assign(paste0(condition, '_leaf_folder'),
           paste0(
             get(sprintf('%s_folders',condition)),
             sprintf('trim_%s/', get(sprintf('%s_Taxon_trim',condition))),
             sprintf('pval_%s/', get(sprintf('%s_SparCC_pval',condition))),
             sprintf('dist_%s_cor_%s/', get(sprintf('%s_Tax_adjacency',condition))
                     , get(sprintf('%s_SparCC_cor',condition))),
             sprintf('fstat_%s/', get(sprintf('%s_Pair_fstat',condition)))
           )
    )
    
    assign(paste0(condition, '_models'),
           list(
             triple=paste0(
               get(sprintf('%s_leaf_folder',condition)),
               sprintf('triplet_models_%s.txt',get(sprintf('%s_Pair_fstat',condition)))
             ),
             double=paste0(
               get(sprintf('%s_leaf_folder',condition)),
               sprintf('pair_models_%s.txt',get(sprintf('%s_Pair_fstat',condition)))
             )
           )
    )
    
    assign(
      paste0(condition, '_triple'),
      read.table(file = get(paste0(condition,'_models'))['triple'][[1]], sep = '\t', header = T, stringsAsFactors = F)
    )
    
    assign(
      paste0(condition, '_double'),
      read.table(file = get(paste0(condition,'_models'))['double'][[1]], sep = '\t', header = T, stringsAsFactors = F)
    )
    
    assign(
      paste0(condition,'_sample'),
      read.table(
        file = paste0(
          get(paste0(condition,'_folders')),
          sprintf('trim_%s/',get(sprintf('%s_Taxon_trim',condition))),
          sprintf('%s.txt',test_or_valid)
        ),
        sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
    )
    
    assign(
      paste0(condition,'_train'),
      read.table(
        file = paste0(
          get(paste0(condition,'_folders')),
          sprintf('trim_%s/',get(sprintf('%s_Taxon_trim',condition))),
          'train.txt'
        ),
        sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
    )
    
    assign(
      paste0(condition,'_tax_code'),
      read.table(
        file = paste0(
          get(paste0(condition,'_folders')),
          sprintf('trim_%s/',get(sprintf('%s_Taxon_trim',condition))),
          'tax_code.txt'
        ),
        sep = '\t', header = F, stringsAsFactors = F)
    )
  }
  
  # Replace taxon names in all data frames with common tax IDs
  tax_code = common_tax_code(I_tax_code,H_tax_code)
  
  for (condition in c('H','I')) {
    temp = get(paste0(condition,'_double'))
    for (i in rownames(temp)){
      temp[i, 'Response'] = tax_code[which(tax_code['Taxon'] == temp[i,'Response']),'Common']
      temp[i, 'Predictor'] = tax_code[which(tax_code['Taxon'] == temp[i,'Predictor']),'Common']
      assign(paste0(condition,'_double'), temp)
    }
    
    temp = get(paste0(condition,'_triple'))
    for (i in rownames(temp)){
      temp[i, 'Response'] = tax_code[which(tax_code['Taxon'] == temp[i,'Response']),'Common']
      temp[i, 'Predictor1'] = tax_code[which(tax_code['Taxon'] == temp[i,'Predictor1']),'Common']
      temp[i, 'Predictor2'] = tax_code[which(tax_code['Taxon'] == temp[i,'Predictor2']),'Common']
      assign(paste0(condition,'_triple'), temp)
    }
    
    map <- setNames(tax_code[[sprintf('%s_ID',condition)]], tax_code[['Common']])
    
    # temp = get(paste0(condition,'_train'))
    # temp2 = rownames(temp)
    # rownames(temp) = map[unlist(temp2)]
    # assign(paste0(condition,'_train'), temp)
    # 
    # Dafuq is rong w dis paragraph?
    # It errors with I_sample
    # and works improperly w H.
    # Y?
    
    
    temp = get(paste0(condition,'_train'))
    new = c()
    for (i in rownames(temp)){
      new = c(new, tax_code[which(tax_code[sprintf('%s_ID',condition)] == i),'Common'])
    }
    rownames(temp) = new
    assign(paste0(condition,'_train'), temp)
    
    
    temp = get(paste0(condition,'_sample'))
    new = c()
    for (i in rownames(temp)){
      new = c(new, tax_code[which(tax_code[sprintf('%s_ID',condition)] == i),'Common'])
    }
    rownames(temp) = new
    assign(paste0(condition,'_sample'), temp)
    
  }
  
  # check if names are placed correctly
  H_mod_taxa = unique(Reduce(union, list(H_double[['Predictor']], H_double[['Response']],
                                  H_triple[['Predictor1']],H_triple[['Predictor2']],H_triple[['Response']])))
  H_tr_taxa = unique(union(rownames(H_train),rownames(H_train)))
  H_mod_not_train = unique(setdiff(H_mod_taxa, H_tr_taxa))
  print(H_mod_not_train)
  
  return(
          list(H_Taxon_trim = H_Taxon_trim,
            H_SparCC_pval = H_SparCC_pval,
            H_SparCC_cor = H_SparCC_cor,
            H_Tax_adjacency = H_Tax_adjacency,
            H_Pair_fstat = H_Pair_fstat,
            I_Taxon_trim = I_Taxon_trim,
            I_SparCC_pval = I_SparCC_pval,
            I_SparCC_cor = I_SparCC_cor,
            I_Tax_adjacency = I_Tax_adjacency,
            I_Pair_fstat = I_Pair_fstat,
            RMSE_sign = RMSE_sign,
            H_double = H_double,
            H_triple = H_triple,
            H_train = H_train,
            H_sample = H_sample,
            I_double = I_double,
            I_triple = I_triple,
            I_train = I_train,
            I_sample = I_sample,
            tax_code = tax_code
            )   
        )
}


# Integrating is pain that is not worth it
# Leavin it here just show I've had this idea
get_pdf_pvals <- function(models,train,test){
  
  calc_err <- function(model, sample){
    
    if (ncol(model) == 5){
      pred = sample[model$Predictor,]
      resp = sample[model$Response,]
      c1 = model$Coef_p
      c2 = model$Intercept
      
      err = (pred*c1 + c2) - resp
    } else {
      pred1 = sample[model$Predictor1,]
      pred2 = sample[model$Predictor2,]
      resp = sample[model$Response,]
      c1 = model$Coef_p1
      c2 = model$Coef_p2
      c3 = model$Intercept
      
      err = (pred1*c1 + pred2*c2 + c3) - resp
    }
    return(err)
  }
  
  approx_pdf <-  function(datum){
    
    bands = 50
    bw_delta = 2
    prev_err = 1
    subdiv = 10000
    at =  1e-4
    rt = 1e-3
    for (i in 1:2000){
      
      kde = density(datum, bw = bands)
      pdf = approxfun(kde$x, kde$y, yleft=0, yright=0)
      total_prob = tryCatch({
        integrate(pdf,-Inf,Inf, subdivisions = subdiv, abs.tol = at, rel.tol = rt)
      }, error = function(e) {
        integrate(pdf,-2*max(abs(datum)),2*max(abs(datum)), subdivisions = subdiv, abs.tol = at, rel.tol = rt)
      })
      err = abs(round(total_prob$value,3) - 1)
      print(round(total_prob$value,3))
      print(err)
      if (err < 0.01) {
        break
      } else {
        if (err > prev_err) {
          bw_delta = -0.8*bw_delta
        }
        bands = bands + bw_delta
        prev_err = err
      }
    }
    print(integrate(pdf,-Inf,Inf, subdivisions = subdiv, abs.tol = at, rel.tol = rt))
    return(kde)
  }
  
  get_pvals <- function(test,model,pdf){
    pvals = list()
    for (ts in colnames(test)){
      err = calc_err(model,test[,ts,drop = F])
      subdiv = 10000
      at =  1e-4
      rt = 1e-3
      pvals[ts] = 1 - integrate(pdf,0,err, subdivisions = subdiv, abs.tol = at, rel.tol = rt)$value
    }
    return(pvals)
  }
  
  pval_table = data.frame(matrix(NA, nrow = nrow(models), ncol = ncol(test) + 1))
  colnames(pval_table) = c('Model',colnames(test))
  
  for (m in rownames(pval_table)){
    errs = c()
    model = models[m,]
    for (tr in colnames(train)){
      errs = c(errs, calc_err(model,train[,tr,drop = F]))
    }
    
    errs_dens = approx_pdf(errs)
    
    pdf = approxfun(errs_dens$x, errs_dens$y, yleft=0, yright=0)
    
    pvals = get_pvals(test,model,pdf)
    
    if (ncol(models) == 5) {
      pval_table[m,1] = paste(models[m,'Response'],models[m,'Predictor'],sep = ',')
      pval_table[m,2:ncol(pval_table)] = pvals
    } else {
      pval_table[m,1] = paste(models[m,'Response'],models[m,'Predictor1'],models[m,'Predictor2'],sep = ',')
      pval_table[m,2:ncol(pval_table)] = pvals
    }
  }
  
  return(pval_table)
}

get_sampling_pval <- function(models,train,test){
  
  calc_err <- function(model, sample){
    
    if (ncol(model) == 5){
      pred = sample[model$Predictor,]
      resp = sample[model$Response,]
      c1 = model$Coef_p
      c2 = model$Intercept
      
      err = abs((pred*c1 + c2) - resp)
    } else {
      pred1 = sample[model$Predictor1,]
      pred2 = sample[model$Predictor2,]
      resp = sample[model$Response,]
      c1 = model$Coef_p1
      c2 = model$Coef_p2
      c3 = model$Intercept
      
      err = abs((pred1*c1 + pred2*c2 + c3) - resp)
    }
    return(err)
  }
  
  get_pvals <- function(test,model, errs){
    pvals = list()
    for (ts in colnames(test)){
      err = calc_err(model,test[,ts,drop = F])
      pvals[ts] = sum(rle(errs >= err)$values)/length(errs)
    }
    return(pvals)
  }
  
  pval_table = data.frame(matrix(NA, nrow = nrow(models), ncol = ncol(test) + 1))
  colnames(pval_table) = c('Model',colnames(test))
  
  for (m in rownames(pval_table)){
    errs = c()
    model = models[m,]
    for (tr in colnames(train)){
      errs = c(errs, calc_err(model,train[,tr,drop = F]))
    }
    
    pvals = get_pvals(test,model, errs)
    
    if (ncol(models) == 5) {
      pval_table[m,1] = paste(models[m,'Response'],models[m,'Predictor'],sep = ',')
      pval_table[m,2:ncol(pval_table)] = pvals
    } else {
      pval_table[m,1] = paste(models[m,'Response'],models[m,'Predictor1'],models[m,'Predictor2'],sep = ',')
      pval_table[m,2:ncol(pval_table)] = pvals
    }
  }
  
  pval_table[,2:ncol(pval_table)] = pval_table[,2:ncol(pval_table)]
  pval_table[nrow(pval_table) + 1,1] = 'Average'
  s = colSums(pval_table[1:(nrow(pval_table)-1),2:ncol(pval_table)]) / (nrow(pval_table) -1)
  pval_table[nrow(pval_table),2:ncol(pval_table)] = s
  
  return(pval_table)
}

pval_to_votes <- function(df, thr){
  df = df[1:(nrow(df)-1),]
  # All cells where pvalue is extremely small
  # are supposed to be a failed test
  # eg: an H model gives pval=0.0001. That means the sample is less likely to be from H group
  df[df < thr] = 0
  df[df >= thr] = 1
  df[nrow(df) + 1,1] = 'Votes'
  df[nrow(df),2:ncol(df)] = round( ( colSums(df[1:(nrow(df)-1),2:ncol(df)]) / (nrow(df) - 1) ) * 100, 2)
  return(df)
}

params = read.csv(params,sep='\t')
#rownames(params)
for (i in 1) {
  temp = params[i,]
  parsed = parse_info(temp, I_folders, H_folders, test_or_valid)
  
}

models = parsed$H_double
train = parsed$H_train
test = parsed$H_sample

h_h = get_sampling_pval(models,train,test)
h_h = pval_to_votes(h_h, 0.2)

models = parsed$I_double
train = parsed$I_train
test = parsed$H_sample

i_h = get_sampling_pval(models,train,test)
i_h = pval_to_votes(i_h, 0.2)

models = parsed$H_double
train = parsed$H_train
test = parsed$I_sample

h_i = get_sampling_pval(models,train,test)
h_i = pval_to_votes(h_i, 0.2)

models = parsed$I_double
train = parsed$I_train
test = parsed$I_sample

i_i = get_sampling_pval(models,train,test)
i_i = pval_to_votes(i_i, 0.2)

h_test = data.frame(matrix(ncol = 20, nrow = 2))
rownames(h_test) = c('I_mod', 'H_mod')
colnames(h_test) = colnames(h_h[2:ncol(h_h)])
h_test['I_mod',] = i_h[nrow(i_h),2:ncol(i_h)]
h_test['H_mod',] = h_h[nrow(h_h),2:ncol(h_h)]
h_test['Result',] = h_test['I_mod',] - h_test['H_mod',]
h_test['Result',][h_test['Result',] >= 0] = 'I'
h_test['Result',][h_test['Result',] < 0] = 'H'

i_test = data.frame(matrix(ncol = 20, nrow = 2))
rownames(i_test) = c('I_mod', 'H_mod')
colnames(i_test) = colnames(i_i[2:ncol(i_i)])
i_test['I_mod',] = i_i[nrow(i_i),2:ncol(i_i)]
i_test['H_mod',] = h_h[nrow(h_i),2:ncol(h_i)]
i_test['Result',] = i_test['I_mod',] - i_test['H_mod',]
i_test['Result',][i_test['Result',] >= 0] = 'I'
i_test['Result',][i_test['Result',] < 0] = 'H'
