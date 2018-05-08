import numpy as np
from sklearn import neighbors, metrics
import pandas as pd
import random
from itertools import combinations

N = {x:[] for x in
    [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    }

N = {x:[] for x in
    [1]
    }

cycles = 1

outcomes = pd.DataFrame(index = None, columns= ['N_nei', 'True_100_models', 'Total_100_models', 'Orgs_in_true_models'])

weights = 'distance'
# weights = 'uniform'

cd_dataset = pd.read_table('../test_mir/simul_cond2(cd).txt', header = 0, index_col=0, sep = '\t')
hc_dataset = pd.read_table('../test_mir/simul_cond1(hc).txt', header = 0, index_col=0, sep = '\t')

hc_zero_count = (hc_dataset == 0).astype(int).sum(axis=1)
hc_zero_count = hc_zero_count[hc_zero_count<(0.5*hc_dataset.shape[1])]
hc_dataset = hc_dataset.loc[list(hc_zero_count.index),]

cd_zero_count = (cd_dataset == 0).astype(int).sum(axis=1)
cd_zero_count = cd_zero_count[cd_zero_count<(0.5*cd_dataset.shape[1])]
cd_dataset = cd_dataset.loc[list(cd_zero_count.index),]

taxa = list(set(cd_dataset.index) & set(hc_dataset.index))
all_pairs = list(combinations(taxa, 2))

tru = []
demat = pd.read_table('../test_mir/dematrix_cd.txt', header = None, index_col = None)
for i in demat.index:
    for j in demat.columns:
        if demat.loc[i,j] > 0:
            tru.append('tax_'+str(i)+'_'+str(j))

for n_neighbors in N:

    print(n_neighbors)

    for c in range(cycles):

        results = pd.DataFrame(index = range(len(all_pairs)), columns= ['Predictor', 'Response', 'Test_sens', 'Test_spec', 'Test_acc'])
        for p in range(len(all_pairs)):

            pair = all_pairs[p]

            pred = pair[0]
            resp = pair[1]

            # points
            cd_pred = cd_dataset.loc[pred,]
            cd_resp = cd_dataset.loc[resp,]
            cd_X = [(cd_pred[i], cd_resp[i]) for i in range(len(cd_resp))]

            cd_X_test_indices = random.sample(xrange(len(cd_X)), int((0.30*len(cd_X))))
            cd_X_train_indices = [i for i in xrange(len(cd_X)) if not(i in cd_X_test_indices)]
            cd_X_test = [cd_X[i] for i in cd_X_test_indices]
            cd_X_train = [cd_X[i] for i in cd_X_train_indices]

            hc_pred = hc_dataset.loc[pred,]
            hc_resp = hc_dataset.loc[resp,]
            hc_X = [(hc_pred[i], hc_resp[i]) for i in range(len(hc_resp))]

            hc_X_test_indices = random.sample(xrange(len(hc_X)), int((0.30*len(hc_X))))
            hc_X_train_indices = [i for i in xrange(len(hc_X)) if not(i in hc_X_test_indices)]
            hc_X_test = [hc_X[i] for i in hc_X_test_indices]
            hc_X_train = [hc_X[i] for i in hc_X_train_indices]

            train_X = np.array(cd_X_train + hc_X_train)
            test_X = np.array(cd_X_test + hc_X_test)

            # labels
            train_y = np.array([1]*len(cd_X_train) + [0]*len(hc_X_train))
            test_y = np.array([1]*len(cd_X_test) + [0]*len(hc_X_test))

            knn = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
            knn.fit(train_X, train_y)

            predicted = knn.predict(test_X)
            conf_matrix = metrics.confusion_matrix(test_y, predicted)

            TP = float(conf_matrix[1, 1])
            TN = float(conf_matrix[0, 0])
            FN = float(conf_matrix[1, 0])
            FP = float(conf_matrix[0, 1])
            sens = TP / (TP+FN)
            spec = TN / (TN+FP)
            acc = (TP+TN)/ (TP+ FN + TN+ FP)

            results.loc[p,] = [pred, resp, sens, spec, acc]

        results.loc[:,['Test_sens', 'Test_spec',  'Test_acc']] = results.loc[:,['Test_sens', 'Test_spec',  'Test_acc']] *100
        results['Test_sens'] = results['Test_sens'].apply(lambda x: round(x,2))
        results['Test_spec'] = results['Test_spec'].apply(lambda x: round(x,2))
        results['Test_acc'] = results['Test_acc'].apply(lambda x: round(x,2))


        results['Real'] = 1
        for i in results.index:
            if (results.loc[i, 'Response'] in tru) and (results.loc[i, 'Predictor'] in tru):
                results.loc[i,'Real'] = 2
            elif not((results.loc[i, 'Response'] in tru) or (results.loc[i, 'Predictor'] in tru)):
                results.loc[i, 'Real'] = 0

        # results.to_csv('mir_knn_n%s.txt' % str(n_neighbors), sep='\t')

        results_100 = results[results['Test_acc'] == 100]
        tru_mods = len(results_100[results_100['Real'] == 2])
        total_mods = results_100.shape[0]
        results_tru = results_100[results_100['Real'] == 2]
        num_tru_orgs = len(set(list(results_tru['Predictor']) + list(results_tru['Response'])))

        outcomes = outcomes.append({'N_nei' : n_neighbors, 'True_100_models':tru_mods, 'Total_100_models':total_mods, 'Orgs_in_true_models' : num_tru_orgs}, ignore_index=True)


        outcomes.to_csv('mir_knn_cross_valid.txt', sep='\t', index=None)
