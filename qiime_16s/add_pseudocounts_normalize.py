import copy
import re
import sys

import numpy as np
import pandas as pd


def add_pseudocounts(counts_path):

    def get_low_high_taxes(tax, all):
        hier = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        tax_lvl = str(tax).split('__')[-2][-1] + '__'

        if tax_lvl == 'k__':
            immediate_higher = None
        else:
            immediate_higher = ','.join(tax.split(',')[:-1])

        if tax_lvl == 's__':
            immediate_lower =None
        else:
            lower_taxa = [x.split(',') for x in all if ((tax in x) and (x != tax))]
            if len(lower_taxa) == 0:
                immediate_lower = None
            else:
                immediate_lower = [','.join(x) for x in lower_taxa if (tax_lvl in x[-2])]
                immediate_lower = list(set(immediate_lower))
        return(immediate_higher, immediate_lower)

    counts = pd.read_table(counts_path, sep='\t', index_col=0, header = 0, engine = 'python')
    counts['Mean'] = counts.mean(axis=1)
    all_taxes = list(counts.index)
    for tax in all_taxes:
        print('Calculating pseudocounts for ' + tax)
        print(str(all_taxes.index(tax)+1) + '/' + str(len(all_taxes)))
        higher, lower = get_low_high_taxes(tax, all_taxes)
        if higher is None:
            uptax = counts.ix[tax,:-1]
        else:
            uptax = counts.ix[higher,:-1]
        if not(lower is None):
            n_subtax = len(lower) + 1
        else:
            n_subtax = 1
        counts.ix[tax,:-1] += (( 0.01 * uptax ) / n_subtax) + 1
    counts = counts.ix[:,:-1].round(decimals = 0)
    return(counts)

def main():
    counts_path = sys.argv[1]
    out = sys.argv[2]
    n = sys.argv[3]
    # counts_path = 'combined.txt'
    # out = 'combined_normalized.txt'
    added = add_pseudocounts(counts_path)
    if n == 'div':
        normalized = added/added.loc['k__Bacteria']
    if n == 'ln':
        normalized = added.apply(np.log)
    if n == 'none':
        normalized = added
    normalized.to_csv(out, sep='\t', header=None)

if __name__ == "__main__":
    main()