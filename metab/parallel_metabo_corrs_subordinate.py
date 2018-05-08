import sys
import pandas as pd


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

print(1)
group_pairs = pd.read_table(sys.argv[1], header = None, index_col= None)

print(group_pairs.head())

pairs = get_significantly_cord_taxa('./13679_hc_cor.out', './13679_hc_pvals.txt', './13679_hc_tax_code.txt')
consum_table = pd.read_table('./consum.txt', sep='\t', header = 0, index_col = 0)
prod_table = pd.read_table('./prod.txt', sep='\t', header = 0, index_col = 0)
consum_table = consum_table[consum_table['Total'] != 0]
prod_table = prod_table[prod_table['Total'] != 0]

print(2)

pairs['pos'] = [(x[0].replace(',',';'), x[1].replace(',',';')) for x in pairs['pos'] if ((x[0].replace(',',';') in prod_table.index) and (x[1].replace(',',';') in prod_table.index))]
pairs['non'] = [(x[0].replace(',',';'), x[1].replace(',',';')) for x in pairs['non'] if ((x[0].replace(',',';') in prod_table.index) and (x[1].replace(',',';') in prod_table.index))]

print(3)

for i in group_pairs.index:

    gg = [
        group_pairs.iloc[i,0],
        group_pairs.iloc[i, 1]
    ]


    raw_numbers_pos = pd.DataFrame(data=None, columns = [gg[0], gg[1]], index = range(2*len(pairs['pos'])))
    raw_numbers_non_sig = pd.DataFrame(data=None, columns = [gg[0], gg[1]], index = range(2*len(pairs['non'])))

    if (sum(consum_table.loc[:, gg[0]]) == 0) or (sum(prod_table.loc[:, gg[1]]) == 0):
        continue

    for p in range(len(pairs['pos'])):
        raw_numbers_pos.loc[2 * p, gg[0]] = consum_table.loc[pairs['pos'][p][0], gg[0]]
        raw_numbers_pos.loc[2 * p, gg[1]] = prod_table.loc[pairs['pos'][p][1], gg[1]]
        raw_numbers_pos.loc[1 + 2 * p, gg[0]] = consum_table.loc[pairs['pos'][p][1], gg[0]]
        raw_numbers_pos.loc[1 + 2 * p, gg[1]] = prod_table.loc[pairs['pos'][p][0], gg[1]]
    for p in range(len(pairs['non'])):
        raw_numbers_non_sig.loc[2 * p, gg[0]] = consum_table.loc[pairs['non'][p][0], gg[0]]
        raw_numbers_non_sig.loc[2 * p, gg[1]] = prod_table.loc[pairs['non'][p][1], gg[1]]
        raw_numbers_non_sig.loc[1 + 2 * p, gg[0]] = consum_table.loc[pairs['non'][p][1], gg[0]]
        raw_numbers_non_sig.loc[1 + 2 * p, gg[1]] = prod_table.loc[pairs['non'][p][0], gg[1]]
        pass

    raw_numbers_pos = raw_numbers_pos.loc[(raw_numbers_pos != 0).any(1)].astype(float)
    raw_numbers_non_sig = raw_numbers_non_sig.loc[(raw_numbers_non_sig != 0).any(1)].astype(float)

    cor_pos = raw_numbers_pos.corr(method='pearson').fillna(0).iloc[1,0]
    cor_non = raw_numbers_non_sig.corr(method='pearson').fillna(0).iloc[1,0]

    boot_cors = []
    non_sig = False
    if raw_numbers_pos.shape[1] > 100:
        N = 100
    else:
        N = raw_numbers_pos.shape[1]

    failed = 0
    for b in range(100):
        temp = pd.concat([raw_numbers_pos.iloc[:, 0].sample(n=N, replace=True, axis=0).reset_index(drop=True),
                          raw_numbers_pos.iloc[:, 1].sample(n=100, replace=True, axis=0).reset_index(drop=True)],
                         axis=1)
        # boot_cors.append(temp.corr(method='pearson').fillna(0).iloc[1,0])
        if abs(temp.corr(method='pearson').fillna(0).iloc[1, 0]) > abs(cor_pos):
            failed += 1
            if failed >= 0.05*N:
                print(gg)
                print('Non-signinficant')
                non_sig = True
                break

    if not(non_sig):

        if abs(cor_pos) > abs(cor_non):
            print(gg)
            print(cor_pos)
            print(cor_non)
