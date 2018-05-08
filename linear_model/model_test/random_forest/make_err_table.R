library(plyr)

# args = c('./inner_hc/trim_0.5/','./inner_cd/trim_0.5/')

args = commandArgs(trailingOnly=TRUE)
print(args)

p_heal_folder = args[1]
p_ibd_folder = args[2]

tax_code_h = read.table(paste0(p_heal_folder,'tax_code.txt'), row.names = 1, stringsAsFactors = F)
tax_code_i = read.table(paste0(p_ibd_folder,'tax_code.txt'), row.names = 1, stringsAsFactors = F)
test_set_h = read.table(paste0(p_heal_folder,'test.txt'), header = 1, stringsAsFactors = F)
test_set_i = read.table(paste0(p_ibd_folder,'test.txt'), header = 1, stringsAsFactors = F)
train_set_h = read.table(paste0(p_heal_folder,'train.txt'), header = 1, stringsAsFactors = F)
train_set_i = read.table(paste0(p_ibd_folder,'train.txt'), header = 1, stringsAsFactors = F)
pairs_h = read.table(paste0(p_heal_folder,'all_models.txt'), header = 1, stringsAsFactors = F)
pairs_i = read.table(paste0(p_ibd_folder,'all_models.txt'), header = 1, stringsAsFactors = F)
triplet_h = read.table(paste0(p_heal_folder,'all_triplet_models.txt'), header = 1, stringsAsFactors = F)
triplet_i = read.table(paste0(p_ibd_folder,'all_triplet_models.txt'), header = 1, stringsAsFactors = F)

columns = c('Predictor', 'Response', 'Coef_p', 'Intercept')
pairs_h = pairs_h[,columns]
pairs_i = pairs_i[,columns]

pairs_h['Condition'] = 'Healthy'
pairs_i['Condition'] = 'IBD'
all_pair = rbind(pairs_h, pairs_i)
rownames(all_pair) = paste0('pm_', rownames(all_pair), '_', all_pair[,'Condition'])

triplet_h['Condition'] = 'Healthy'
triplet_i['Condition'] = 'IBD'
all_trip = rbind(triplet_h, triplet_i)
rownames(all_trip) = paste0('tm_', rownames(all_trip), '_', all_trip[,'Condition'])

rownames(test_set_h) = tax_code_h[rownames(test_set_h),]
rownames(test_set_i) = tax_code_i[rownames(test_set_i),]
rownames(train_set_h) = tax_code_h[rownames(train_set_h),]
rownames(train_set_i) = tax_code_i[rownames(train_set_i),]
train_set_h['Condition',] = 'Healthy'
train_set_i['Condition',] = 'IBD'
test_set_h['Cols',] = as.character(colnames(test_set_h))
test_set_i['Cols',] = as.character(colnames(test_set_i))
train_set_h['Cols',] = as.character(colnames(train_set_h))
train_set_i['Cols',] = as.character(colnames(train_set_i))

train_all = data.frame(t(rbind.fill(data.frame(t(train_set_h)), data.frame(t(train_set_i)))),stringsAsFactors = F)
colnames(train_all) = as.character(train_all['Cols',])
train_all['Cols'] <- NULL
train_all[is.na(train_all)] = 0
rownames(train_all) = gsub('\\.', ',', x = rownames(train_all))

errs = data.frame(matrix(NA, ncol = ncol(train_all), nrow = nrow(all_pair) + nrow(all_trip) + 1))
rownames(errs) = c(rownames(all_pair), rownames(all_trip), 'Condition')
colnames(errs) = colnames(train_all)
errs['Condition',] = as.character(train_all['Condition',])

get_err_p <- function(model, sample){
  coef =as.numeric(model[['Coef_p']])
  intercept = as.numeric(model[['Intercept']])
  predictor = as.numeric(sample[model[['Predictor']],])
  response = as.numeric(sample[model[['Response']],])
  prediction = coef*predictor + intercept
  err = prediction - response
  return(err)
}

get_err_t <- function(model, sample){
  coef1 =as.numeric(model[['Coef_p1']])
  coef2 =as.numeric(model[['Coef_p2']])
  intercept = as.numeric(model[['Intercept']])
  predictor1 = as.numeric(sample[model[['Predictor1']],])
  predictor2 = as.numeric(sample[model[['Predictor2']],])
  response = as.numeric(sample[model[['Response']],])
  prediction = coef1*predictor1 + coef2*predictor2 + intercept
  err = prediction - response
  return(err)
}


# It's gonna take a while.
# Most likely a whole day.
# So feel fre to parallelize this if you don't have the time

for (mod in rownames(errs)) {
  print(mod)
  if (! is.na(errs[mod,1])) {
    next
  }
  print(which(rownames(errs) == mod))
  if (mod == 'Condition') {
    next
  }
  for (sam in colnames(errs)) {
    if (substr(mod, 1,2) == 'pm') {
      errs[mod, sam] = get_err_p(all_pair[mod,],train_all[sam])
    }
    if (substr(mod, 1,2) == 'tm') {
      errs[mod, sam] = get_err_t(all_trip[mod,],train_all[sam])
    }
  }
}

errs$Rows = rownames(errs)

errs_complete = errs[complete.cases(errs), ]

if ('Rows' %in% colnames(errs_complete)){
  errs_complete$Rows = NULL
}

# write.table makes a left-shifted header
# write.csv does not. However, it can't set separator in large tables.
write.table(errs_complete, file = 'all_errs.txt', quote = F, sep = '\t', col.names=NA)
write.table(all_pair, file = 'all_pair_models.txt', quote = F, sep = '\t')
write.table(all_trip, file = 'all_triple_models.txt', quote = F, sep = '\t')
