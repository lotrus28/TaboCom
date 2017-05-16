import sys
import pandas as pd
from io import StringIO

def parse_collapsed(path_counts):

    print('Reading file: ' + path_counts)
    with open(path_counts) as f:
        input = StringIO(f.read().replace('", ""', '').replace('"', '').replace(', ', ',').replace('\0',''))
        counts = pd.read_table(input, sep='\t', index_col=2, header=None, engine='python')

    counts = counts[1]
    counts = counts.groupby(counts.index).sum()
    return(counts)

def add_higher_taxa(counts):
    hier = ['s__', 'g__', 'f__', 'o__', 'c__', 'p__', 'k__']
    new_counts = counts.copy()
    for lvl in hier:
        this_lvl_taxa = [x for x in counts.index
                         if (lvl in x.split(sep=',')[-1])]
        for tax in this_lvl_taxa:
            tax_terms = tax.split(sep=',')
            higher_taxa = [','.join(tax_terms[:i + 1]) for i in range(len(tax_terms) - 1)]
            for higher_tax in higher_taxa:
                if higher_tax in new_counts.index:
                    new_counts[higher_tax] += counts[tax]
                else:
                    new_counts[higher_tax] = counts[tax]
    return(new_counts)

def reformat_count_table(df):
    if float(pd.__version__[:4]) >= 0.17:
        df = df.sort_values(ascending=False)
    else:
        df = df.sort(ascending=False)

    output = pd.DataFrame(data=0, index=range(len(df)), columns=['Counts', 'Taxon'])
    output['Counts'] = df.values
    output['Taxon'] = df.index
    return(output)

def main():
    path_counts = sys.argv[1]
    path_out = sys.argv[2]
    counts = parse_collapsed(path_counts)
    higher_added = add_higher_taxa(counts)
    out = reformat_count_table(higher_added)
    path_counts = path_counts.split('/')[-1]
    out.to_csv(path_out + '/' + path_counts.replace('.txt', '_added_higher.txt'), header=None, index= None, sep='\t')

if __name__ == "__main__":
    main()