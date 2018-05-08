import pandas as pd
from itertools import combinations
import subprocess
import time

def get_significantly_cord_taxa(p_cor, p_pval, p_taxa):

    hc_cor = pd.read_table(p_cor, sep='\t', header=0, index_col=0)
    hc_pval = pd.read_table(p_pval, sep='\t', header=0, index_col=0)
    hc_taxa = pd.read_table(p_taxa, sep='\t', header=None, index_col=0)

    hc_cor.index = [hc_taxa.loc[x, 1] for x in hc_cor.index]
    hc_cor.columns = [hc_taxa.loc[x, 1] for x in hc_cor.columns]
    hc_pval.index = [hc_taxa.loc[x, 1] for x in hc_pval.index]
    hc_pval.columns = [hc_taxa.loc[x, 1] for x in hc_pval.columns]

    pos_sig_pairs = []
    neg_sig_pairs = []
    non_sig_pairs = []
    for i in hc_cor.index[:-1]:
        for j in hc_cor.columns[list(hc_cor.index).index(i) + 1:]:

            if hc_pval.loc[i, j] < 0.001:

                if hc_cor.loc[i, j] >= 0.3:
                    pos_sig_pairs.append((i, j))
                elif hc_cor.loc[i, j] <= (-0.3):
                    print('11111111111111111111')
                    neg_sig_pairs.append((i, j))
                else:
                    non_sig_pairs.append((i, j))
            else:
                non_sig_pairs.append((i, j))

    return({'pos':pos_sig_pairs, 'neg':neg_sig_pairs, 'non':non_sig_pairs})

nodes = [
"node15.net0.pyxis.ripcm.com",
"node5.net0.pyxis.ripcm.com",
"node8.net0.pyxis.ripcm.com",
"node10.net0.pyxis.ripcm.com",
"node12.net0.pyxis.ripcm.com",
"node14.net0.pyxis.ripcm.com"
]

n_streams = len(nodes)


# Check if consum/prod are correlated
consum_table = pd.read_table('./consum.txt', sep='\t', header = 0, index_col = 0)
prod_table = pd.read_table('./prod.txt', sep='\t', header = 0, index_col = 0)
consum_table = consum_table[consum_table['Total'] != 0]
prod_table = prod_table[prod_table['Total'] != 0]

pairs = get_significantly_cord_taxa('./13679_hc_cor.out', './13679_hc_pvals.txt', './13679_hc_tax_code.txt')

pairs['pos'] = [(x[0].replace(',',';'), x[1].replace(',',';')) for x in pairs['pos'] if ((x[0].replace(',',';') in prod_table.index) and (x[1].replace(',',';') in prod_table.index))]
pairs['non'] = [(x[0].replace(',',';'), x[1].replace(',',';')) for x in pairs['non'] if ((x[0].replace(',',';') in prod_table.index) and (x[1].replace(',',';') in prod_table.index))]

groups = pd.read_table('./consum.txt', sep='\t', header = 0, index_col = 0)
group_names = list(groups.columns)[:-1]
group_pairs = list(combinations(group_names, 2))

stream_group_pairs = [group_pairs[i:i +len(group_pairs)/ n_streams] for i in range(0, len(group_pairs), len(group_pairs)/n_streams)]
if len(stream_group_pairs[-1]) != len(stream_group_pairs[0]):
    stream_group_pairs[-2] = stream_group_pairs[-2] + stream_group_pairs[-1]
    stream_group_pairs = stream_group_pairs[:-1]
stream_group_pairs = {'stream_%s.txt'%i : stream_group_pairs[i] for i in range(len(stream_group_pairs))}

for s in stream_group_pairs:
    time.sleep(1)
    with open(s, 'a') as f:
        for pair in stream_group_pairs[s]:
            f.write("%s\t%s\n" % (pair[0], pair[1]))

    command = 'echo "python calc_metab_cors.py %s >> %s.log" | qsub -pe make 3 -N stream%s -cwd -l hostname=%s' % (s,s, int(s.split('_')[1][0]), nodes[int(s.split('_')[1][0])])
    print(command)
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

