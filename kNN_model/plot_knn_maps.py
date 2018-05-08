from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from sklearn import neighbors, metrics
import numpy as np
import pandas as pd
import random as random
import time

# plotting
# step size in the mesh
h = 1
n_neighbors = 5
weights = 'distance'

# Create color maps
cmap_light = ListedColormap(['#FFAAAA', '#AAFFAA'])
cmap_bold = ListedColormap(['#FF0000', '#00FF00'])


cd_dataset = pd.read_table('./prjeb/13679_otu_ibd.txt', header = 0, index_col=0, sep = '\t')
hc_dataset = pd.read_table('./prjeb/13679_otu_hc.txt', header = 0, index_col=0, sep = '\t')

hc_zero_count = (hc_dataset == 0).astype(int).sum(axis=1)
hc_zero_count = hc_zero_count[hc_zero_count<(0.5*hc_dataset.shape[1])]
hc_dataset = hc_dataset.loc[list(hc_zero_count.index),]

cd_zero_count = (cd_dataset == 0).astype(int).sum(axis=1)
cd_zero_count = cd_zero_count[cd_zero_count<(0.5*cd_dataset.shape[1])]
cd_dataset = cd_dataset.loc[list(cd_zero_count.index),]

common_taxa = list(set(cd_dataset.index) & set(hc_dataset.index))
cd_dataset = cd_dataset.loc[common_taxa,]
hc_dataset = hc_dataset.loc[common_taxa,]

# all_pairs = list(combinations(common_taxa, 2))

all_pairs = []
select_class = pd.read_table('./prjeb/prjeb_knn_2_more67.txt', header = 0, index_col = 0, sep = '\t')
for i in select_class.index:
    all_pairs.append((select_class.loc[i,'Predictor'],select_class.loc[i,'Response']))

results = pd.DataFrame(index = range(len(all_pairs)), columns= ['Predictor', 'Response', 'Test_sens', 'Test_spec', 'Test_acc', 'Max_Individual'])

for p in range(len(all_pairs))[2:]:
    pair = all_pairs[p]

    pred = pair[0]
    resp = pair[1]

    # points
    cd_pred = cd_dataset.loc[pred,]
    cd_resp = cd_dataset.loc[resp,]
    cd_X = [(cd_pred[i], cd_resp[i]) for i in range(len(cd_resp))]

    cd_X_test_indices = random.sample(xrange(len(cd_X)), int((0.30 * len(cd_X))))
    cd_X_train_indices = [i for i in xrange(len(cd_X)) if not (i in cd_X_test_indices)]
    cd_X_test = [cd_X[i] for i in cd_X_test_indices]
    cd_X_train = [cd_X[i] for i in cd_X_train_indices]

    hc_pred = hc_dataset.loc[pred,]
    hc_resp = hc_dataset.loc[resp,]
    hc_X = [(hc_pred[i], hc_resp[i]) for i in range(len(hc_resp))]

    hc_X_test_indices = random.sample(xrange(len(hc_X)), int((0.30 * len(hc_X))))
    hc_X_train_indices = [i for i in xrange(len(hc_X)) if not (i in hc_X_test_indices)]
    hc_X_test = [hc_X[i] for i in hc_X_test_indices]
    hc_X_train = [hc_X[i] for i in hc_X_train_indices]

    train_X = np.array(cd_X_train + hc_X_train)
    test_X = np.array(cd_X_test + hc_X_test)

    # labels
    train_y = np.array([1] * len(cd_X_train) + [0] * len(hc_X_train))
    test_y = np.array([1] * len(cd_X_test) + [0] * len(hc_X_test))

    # we create an instance of Neighbours Classifier and fit the data.

    knn = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
    knn.fit(train_X, train_y)

    predicted = knn.predict(test_X)
    conf_matrix = metrics.confusion_matrix(test_y, predicted)

    TP = float(conf_matrix[1, 1])
    TN = float(conf_matrix[0, 0])
    FN = float(conf_matrix[1, 0])
    FP = float(conf_matrix[0, 1])
    sens = TP / (TP + FN)
    spec = TN / (TN + FP)
    acc = (TP + TN) / (TP + FN + TN + FP)

    acc = round(acc * 100, 2)

    # X = test_X
    # y = test_y

    # X = train_X
    # y = train_y

    X = np.array(cd_X_train + hc_X_train + cd_X_test + hc_X_test)
    y = np.array([1] * len(cd_X_train) + [0] * len(hc_X_train) + [1] * len(cd_X_test) + [0] * len(hc_X_test))

    # Plot the decision boundary. For that, we will assign a color to each
    # point in the mesh [x_min, x_max]x[y_min, y_max].
    # x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    # y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    x_min, x_max = X[:, 0].min() - 1, np.percentile(X[:, 0], 90) + 1
    y_min, y_max = X[:, 1].min() - 1, np.percentile(X[:, 1], 90) + 1



    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))
    Z = knn.predict(np.c_[xx.ravel(), yy.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)

    fig = plt.figure()
    plt.pcolormesh(xx, yy, Z, cmap=cmap_light)

    # Plot also the training pointss
    plt.scatter(X[:, 0], X[:, 1], c=y, cmap=cmap_bold,
                edgecolor='k', s=20)
    plt.xlim(xx.min(), xx.max())
    plt.ylim(yy.min(), yy.max())
    plt.title("2-Class classification (k = %i, weights = '%s', acc = %s" % (n_neighbors, weights, acc) + "%)"
              )
    plt.xlabel(pred)
    plt.ylabel(resp)
    plt.savefig('./%s.png' % p)
    plt.close(fig)
    time.sleep(1)


###################################################