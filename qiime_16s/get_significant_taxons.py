import sys
from difflib import SequenceMatcher
import pandas as pd

global cor_thr
cor_thr = float(sys.argv[3])
global tax_thr
tax_thr = int(float(sys.argv[4]))

taxons = ['g__', 'f__', 'o__', 'c__', 'p__']
tax_name = taxons[tax_thr]

# in: sig.cors table from R
# out: same table w.removed positive cor-s for close taxons and taxon-subtaxon cor-s
def get_edges(sig_cors, out):
    t = pd.read_table(sig_cors, header=0, index_col=0, engine='python')
    sig_edges = pd.DataFrame(columns=range(3))

    for row in list(t.index)[:-1]:
        i = list(t.index).index(row)
        for col in list(t.columns)[i + 1:]:
            # All strong negative correlations are significant
            if t.loc[row, col] < -(cor_thr):
                temp = pd.DataFrame([row, col, t.loc[row, col]], columns=range(1))
                sig_edges = sig_edges.append(temp.transpose())
            if t.loc[row, col] > cor_thr:
                temp = SequenceMatcher(None, col, row)
                temp = temp.find_longest_match(0, len(col), 0, len(row))
                # If the common taxonomy part is not 1 organism's full taxonomy
                if not (temp[2] == min(len(col), len(row))):
                    lcs = col[temp[0]:temp[0] + temp[2]]
                    # Strong positive cor-s are significant
                    # only if they occur at order level or higher
                    if not (tax_name in lcs):
                        temp = pd.DataFrame([row, col, t.loc[row, col]], columns=range(1))
                        sig_edges = sig_edges.append(temp.transpose())

    sig_edges = sig_edges.sort_values(2,0,ascending = False)

    outfile = out + '/edges_' + str(cor_thr) + '_' + str(float(tax_thr)) + '.txt'
    sig_edges.to_csv(outfile, sep='\t', header=False, index=False)

    if sig_edges.shape[0] == 0:
        with open(outfile,'w') as f:
            f.write('No edges based on pair correlations\n')
            print('Empty file written')
    return ()

print('get_significant_taxons.py: ' + str(cor_thr) + '__' + str(tax_thr))
get_edges(sys.argv[1], sys.argv[2])