'apteka.mbdb@gmail.com'
import sys
import pandas

def combine_otu_tables(path_to_otu):
    t = pandas.read_table(path_to_otu, header=None)
    # filter out rare bacteria
    t = t[(t == 0).astype(int).sum(axis=1) < t.shape[1]/2]
    t.columns = ['pat' + str(x) for x in t.columns]
    t.index = ['tax' + str(x+1) for x in t.index]
    tax_code = t.iloc[:,0]
    t = t.drop('pat0', 1)
    t.iloc[:, :(t.shape[1]-20)].to_csv('train.txt', sep = '\t')
    t.iloc[:, (t.shape[1] - 20):].to_csv('test.txt', sep='\t')
    tax_code.to_csv('tax_code.txt', sep = '\t')
    return()

combine_otu_tables(sys.argv[1])