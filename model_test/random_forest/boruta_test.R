library(Boruta)
library(randomForest)

errs = read.table('errs_to_boruta.txt', header = 1, row.names = 1)
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

feat = Boruta(x = train_errs, y = conds[[1]], maxRuns = 200, doTrace = 2, 
              holdHistory = F, getImp = getImpLegacyRfZ, num.trees = 10000)
conf = names(feat$finalDecision[feat$finalDecision == 'Confirmed'])
tent = names(feat$finalDecision[feat$finalDecision == 'Tentative'])
# rej = names(feat$finalDecision[feat$finalDecision == 'Rejected'])
mRuns = 100
for (i in 1:mRuns) {
  if (length(tent) == 0) {
    break
  }
  feat = Boruta(x = train_errs[,c(conf, tent)], y = conds[[1]], maxRuns = 200, doTrace = 2, 
                holdHistory = F, getImp = getImpLegacyRfZ, num.trees = 10000)
  conf = unique(c(conf,names(feat$finalDecision[feat$finalDecision == 'Confirmed'])))
  tent = names(feat$finalDecision[feat$finalDecision == 'Tentative'])
}


# conds = train_errs['Condition']
# train_errs$Condition = NULL
# classifier <-  randomForest(x = train_errs, y = conds, ntree = 3*length(conf))

train_errs['Condition'] = conds
train_errs = train_errs[c(conf,'Condition')]
classifier <-  randomForest(Condition ~ ., data = train_errs, ntree = 3*length(conf))


tesr_errs_cond = test_errs['Condition']
test_errs$Condition = NULL
preds = predict(classifier, test_errs[conf])

result=cbind(tesr_errs_cond, preds)
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

sens = 100 * nrow(result[result['Quality'] == 'TP']) / (nrow(result[result['Quality'] == 'TP']) + nrow(result[result['Quality'] == 'FN']))
spec = 100 * nrow(result[result['Quality'] == 'TN']) / (nrow(result[result['Quality'] == 'TN']) + nrow(result[result['Quality'] == 'FP']))

print( paste0('Sensitivity: ', round(sens,2)) )
print( paste0('Specificity: ', round(spec,2)) )

write(conf,'conf_models.txt')
write(tent,'tent_models.txt')
save(classifier, file = 'randomForest.Rdata')