import pandas as pd
import sys

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
                    neg_sig_pairs.append((i, j))
                else:
                    non_sig_pairs.append((i, j))
            else:
                non_sig_pairs.append((i, j))

    return({'pos':pos_sig_pairs, 'neg':neg_sig_pairs, 'non':non_sig_pairs})


prefix = sys.argv[1]

with open('./stream0.txt') as f:
    group_pairs = f.read().splitlines()

pairs = get_significantly_cord_taxa('./%s_cor.out'%prefix, './%s_pvals.txt'%prefix, './%s_tax_code.txt'%prefix)
consum_table = pd.read_table('./consum.txt', sep='\t', header = 0, index_col = 0)
prod_table = pd.read_table('./prod.txt', sep='\t', header = 0, index_col = 0)
consum_table = consum_table[consum_table['Total'] != 0]
prod_table = prod_table[prod_table['Total'] != 0]

pairs['pos'] = [(x[0].replace(',',';'), x[1].replace(',',';')) for x in pairs['pos'] if ((x[0].replace(',',';') in prod_table.index) and (x[1].replace(',',';') in prod_table.index))]
pairs['non'] = [(x[0].replace(',',';'), x[1].replace(',',';')) for x in pairs['non'] if ((x[0].replace(',',';') in prod_table.index) and (x[1].replace(',',';') in prod_table.index))]


counter = 0
for i in group_pairs:

    counter += 1
    if counter % 5 == 1:
        print('%s out of %s' % (counter, len(group_pairs)))

    gg = [i,i]

    raw_numbers_pos = pd.DataFrame(data=None, columns=[0, 1], index=range(2 * len(pairs['pos'])))
    raw_numbers_non_sig = pd.DataFrame(data=None, columns=[0, 1], index=range(2 * len(pairs['non'])))

    if (sum(consum_table.loc[:, gg[0]]) == 0) or (sum(prod_table.loc[:, gg[1]]) == 0):
        continue

    for p in range(len(pairs['pos'])):
        raw_numbers_pos.iloc[2 * p, 0] = consum_table.loc[pairs['pos'][p][0], gg[0]]
        raw_numbers_pos.iloc[2 * p, 1] = prod_table.loc[pairs['pos'][p][1], gg[1]]
        raw_numbers_pos.iloc[1 + 2 * p, 0] = consum_table.loc[pairs['pos'][p][1], gg[0]]
        raw_numbers_pos.iloc[1 + 2 * p, 1] = prod_table.loc[pairs['pos'][p][0], gg[1]]
    for p in range(len(pairs['non'])):
        raw_numbers_non_sig.iloc[2 * p, 0] = consum_table.loc[pairs['non'][p][0], gg[0]]
        raw_numbers_non_sig.iloc[2 * p, 1] = prod_table.loc[pairs['non'][p][1], gg[1]]
        raw_numbers_non_sig.iloc[1 + 2 * p, 0] = consum_table.loc[pairs['non'][p][1], gg[0]]
        raw_numbers_non_sig.iloc[1 + 2 * p, 1] = prod_table.loc[pairs['non'][p][0], gg[1]]

    raw_numbers_pos = raw_numbers_pos.loc[(raw_numbers_pos != 0).any(1)].astype(float)
    raw_numbers_non_sig = raw_numbers_non_sig.loc[(raw_numbers_non_sig != 0).any(1)].astype(float)

    if raw_numbers_pos.shape[0] > 100:
        N = 100
    else:
        N = raw_numbers_pos.shape[0]

    if N < 4:
        continue

    cor_pos = raw_numbers_pos.corr(method='pearson').fillna(0).iloc[1, 0]
    cor_non = raw_numbers_non_sig.corr(method='pearson').fillna(0).iloc[1, 0]

    if abs(cor_pos) > abs(cor_non):

        boot_cors = []
        non_sig = False

        failed = 0
        for b in range(1000):
            temp = pd.concat([raw_numbers_pos.iloc[:, 0].sample(n=N, replace=True, axis=0).reset_index(drop=True),
                              raw_numbers_pos.iloc[:, 1].sample(n=N, replace=True, axis=0).reset_index(drop=True)],
                             axis=1)
            # boot_cors.append(temp.corr(method='pearson').fillna(0).iloc[1,0])
            if abs(temp.corr(method='pearson').fillna(0).iloc[1, 0]) > abs(cor_pos):
                failed += 1

        pval = round(float(failed)/10.,2)

        with open('%s_metab_pv.out'%prefix, 'a') as f:
            f.write("%s\t%s\t%s\t%s'\n'" % (gg[1], cor_pos, cor_non, pval))