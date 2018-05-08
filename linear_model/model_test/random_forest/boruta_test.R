library(Boruta)
library(randomForest)

args = commandArgs(trailingOnly=TRUE)
print(args)

#args = '/data5/bio/runs-galkin/PRJEB13679/log_errs/all_errs.txt'

p_err_tabe = args[1]
# p_err_table = 'errs_to_boruta.txt'
errs = read.table(p_err_tabe, header = 1, row.names = 1, stringsAsFactors = F)
# errs = as.data.frame(t(errs))

if ('Condition' %in% rownames(errs)) {
	errs = as.data.frame(t(errs))
}

ibd_errs = errs[errs['Condition'] == 'IBD',]
heal_errs = errs[errs['Condition'] == 'Healthy',]

test_errs = merge(x = ibd_errs[1:(nrow(ibd_errs)*0.15),], y = heal_errs[1:(nrow(heal_errs)*0.15),], all = TRUE)
rownames(test_errs) = c(rownames(ibd_errs[1:(nrow(ibd_errs)*0.15),]), rownames(heal_errs[1:(nrow(heal_errs)*0.15),]))
train_errs = merge(x = ibd_errs[((nrow(ibd_errs)*0.15)+1):nrow(ibd_errs),], y = heal_errs[((nrow(heal_errs)*0.15)+1):nrow(heal_errs),], all = TRUE)
rownames(train_errs) = c(rownames(ibd_errs[((nrow(ibd_errs)*0.15)+1):nrow(ibd_errs),]), rownames(heal_errs[((nrow(heal_errs)*0.15)+1):nrow(heal_errs),]))

rm(heal_errs)
rm(ibd_errs)
conds = train_errs['Condition']
train_errs$Condition = NULL

for (i in colnames(test_errs)){
  
  if (i != 'Condition'){
    test_errs[,i] = as.numeric(as.character(test_errs[,i]))
  }
}

for (i in colnames(train_errs)){
  train_errs[,i] = as.numeric(as.character(train_errs[,i]))
}

feat = Boruta(x = train_errs, y = factor(conds[[1]]), maxRuns = 4000, doTrace = 2, 
              holdHistory = F, getImp = getImpLegacyRfZ, num.trees = 10000)
conf = names(feat$finalDecision[feat$finalDecision == 'Confirmed'])
tent = names(feat$finalDecision[feat$finalDecision == 'Tentative'])
# rej = names(feat$finalDecision[feat$finalDecision == 'Rejected'])

# if (length(c(conf,tent)) != 0) {
#   mRuns = 10
#   for (i in 1:mRuns) {
#     if (length(tent) != 0) {
#       feat = Boruta(x = train_errs[,c(conf, tent)], y = factor(conds[[1]]), maxRuns = 200, doTrace = 2, 
#                     holdHistory = F, getImp = getImpLegacyRfZ, num.trees = 10000)
#       conf = unique(c(conf,names(feat$finalDecision[feat$finalDecision == 'Confirmed'])))
#       tent = names(feat$finalDecision[feat$finalDecision == 'Tentative'])
#     }
#   }
# }

# conds = train_errs['Condition']
# train_errs$Condition = NULL
# classifier <-  randomForest(x = train_errs, y = factor(conds), ntree = 3*length(conf))

train_errs['Condition'] = factor(conds$Condition)
train_errs = train_errs[c(conf,'Condition')]
classifier <-  randomForest(Condition ~ ., data = train_errs, ntree = 3*length(conf))


test_errs_cond = test_errs['Condition']
test_errs$Condition = NULL
test_errs[] <- as.numeric(as.character(as.matrix(test_errs)))
preds = predict(classifier, test_errs[conf])

result=cbind(test_errs_cond, preds)
result$Quality = NA
for (i in rownames(result)){
  if (result[i,'Condition'] == 'Healthy' && result[i,'preds'] == 'Healthy') {
    result[i,'Quality'] = 'TN'
  }
  if (result[i,'Condition'] == 'IBD' && result[i,'preds'] == 'IBD') {
    result[i,'Quality'] = 'TP'
  }
  if (result[i,'Condition'] == 'IBD' && result[i,'preds'] == 'Healthy') {
    result[i,'Quality'] = 'FN'
  }
  if (result[i,'Condition'] == 'Healthy' && result[i,'preds'] == 'IBD') {
    result[i,'Quality'] = 'FP'
  }
}


sens = 100 * nrow(result[result['Quality'] == 'TP',]) / (nrow(result[result['Quality'] == 'TP',]) + nrow(result[result['Quality'] == 'FN',]))
spec = 100 * nrow(result[result['Quality'] == 'TN',]) / (nrow(result[result['Quality'] == 'TN',]) + nrow(result[result['Quality'] == 'FP',]))

print( paste0('Sensitivity: ', round(sens,2)) )
print( paste0('Specificity: ', round(spec,2)) )

write(conf,'conf_models.txt')
write(tent,'tent_models.txt')
save(classifier, file = 'randomForest.Rdata')
