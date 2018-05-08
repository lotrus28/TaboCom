import sys
import pandas

global Taxon_trim
Taxon_trim = float(sys.argv[3])

# in: otu table for many patients
# out: training and testing set otus with reduced taxon names.
def combine_otu_tables(path_to_otu, out):

    t = pandas.read_table(path_to_otu, header=0, engine = 'python')
    # Replce taxon names with IDS
    # t.columns = ['pat' + str(x) for x in t.columns]
    t.index = ['tax' + str(x+1) for x in t.index]
    # Remember IDs
    tax_code = t.iloc[:,0]
    t = t.drop(t.columns[0], 1)

    train_all = t.iloc[:, :(t.shape[1]-40)]

    # Filter out rare bacteria
    # Threshold is set in teach_models.sh
    # Possible values are: 0.1 - 1.0
    train_filt = train_all[(train_all == 0).astype(int).sum(axis=1) < (train_all.shape[1]*Taxon_trim)]

    # Write down training sets (both filtered and not)
    train_filt.to_csv(out + '/train.txt', sep='\t')
    train_all.to_csv(out + '/train_all_taxons.txt', sep = '\t')
    # Write down non-filtered testing set
    t.iloc[:, (t.shape[1] - 40):(t.shape[1] - 20)].to_csv(out +'/test.txt', sep='\t')
    # Write down non-filtered validation set
    t.iloc[:, (t.shape[1] - 20):].to_csv(out + '/valid.txt', sep='\t')
    tax_code.to_csv(out + '/tax_code.txt', sep = '\t')
    return()

print('sep_sets.py: ' + str(Taxon_trim))
combine_otu_tables(sys.argv[1], sys.argv[2])