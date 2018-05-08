import pandas as pd
import sys
def scramble_train_sets(path_to_otus,out = '.'):
    abun_tables = []
    total_index = []
    for i in path_to_otus:
        abun_tables.append(
            pd.read_table(i, header=0, index_col=0, engine='python')
        )
        total_index += list(abun_tables[-1].index)
    for otu_tab in abun_tables:
        for tax in total_index:
            if not tax in otu_tab.index:
                otu_tab.ix[tax,:] = 0
    train_sets = [x.iloc[:, :(x.shape[1] - 40)] for x in abun_tables]
    test_sets = [x.iloc[:, (x.shape[1] - 40):] for x in abun_tables]
    combined = pd.concat(train_sets, axis=1)
    L = len(train_sets)
    scrambled_train = [None for _ in range(L)]
    for i in range(combined.shape[1]):
        j = i%L
        if scrambled_train[j] is None:
            scrambled_train[j] = combined.iloc[:,i]
            continue
        scrambled_train[j] = pd.concat([scrambled_train[j],combined.iloc[:,i]], axis = 1)
    scr_train__saved_test = [pd.concat([scrambled_train[i],test_sets[i]], axis = 1) for i in range(L)]
    scr_train__saved_test = [x.loc[(x != 0).any(1)] for x in scr_train__saved_test]
    i = 0
    for el in scr_train__saved_test:
        el.to_csv('{}/scrambled_{}.txt'.format(out,str(i)), sep='\t')
        i+=1
    return ()

def main():
    x = sys.argv[1]
    y = sys.argv[2]
    try:
        out = sys.argv[3]
    except:
        out = '.'
    scramble_train_sets([x, y],out)

if __name__ == '__main__':
    main()
