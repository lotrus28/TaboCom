import sys
import re
import pandas as pd

def combine_otu_tables(path_to_files):

    with open(path_to_files) as a:
        filenames = a.read().splitlines()

    separated = {re.search(r'ERR\d+?(?=_)',x).group(0):pd.read_table(x, sep = '\t', index_col = 1, header = None,engine='python')
                 for x in filenames}
    indices = [list(x.index) for x in list(separated.values())]
    all_taxa = sum(indices,[])
    all_taxa = list(set(all_taxa))
    altogether = pd.DataFrame(None, columns = list(separated.keys()), index = all_taxa)

    for pat in separated:
        altogether[pat] = separated[pat][0]
    altogether = altogether.fillna(0)

    altogether['Mean'] = altogether.mean(axis = 1)

    if float(pd.__version__[:4]) >= 0.17:
        altogether = altogether.sort_values('Mean', axis = 0, ascending=False)
    else:
        altogether = altogether.sort('Mean', axis = 0, ascending=False)

    return(altogether.ix[:,:-1])

def main():
    # list_of_files = 'temp2.txt'
    # output = 'combined.txt'
    list_of_files = sys.argv[1]
    output = sys.argv[2]
    combined = combine_otu_tables(list_of_files)
    print('Combining all OTU-tables')
    combined.to_csv(output, sep = '\t')

if __name__ == "__main__":
    main()
