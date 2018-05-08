import pandas as pd
from itertools import combinations

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
# group_pairs = list(combinations(group_names, 2))
group_pairs = [(x,x) for x in group_names]

print(len(group_pairs))

counter = 0
for gg in group_pairs:

    counter += 1
    if counter%10 == 1:
        print('%s out of %s'%(counter , len(group_pairs)))

    raw_numbers_pos = pd.DataFrame(data=None, columns = [0,1], index = range(2*len(pairs['pos'])))
    raw_numbers_non_sig = pd.DataFrame(data=None, columns = [0,1], index = range(2*len(pairs['non'])))

    if (sum(consum_table.loc[:, gg[0]]) == 0) or (sum(prod_table.loc[:, gg[1]]) == 0):
        continue
	
	# Each correlated pair has 2 pairs of (cons.score,prod.score) values
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


    cor_pos = raw_numbers_pos.corr(method='pearson').fillna(0).iloc[1,0]
    cor_non = raw_numbers_non_sig.corr(method='pearson').fillna(0).iloc[1,0]

    if raw_numbers_pos.shape[0] > 100:
        N = 100
    else:
        N = raw_numbers_pos.shape[0]

    if N < 4:
        continue
	
	# Sign.corr-n is supp-d to be stronger than in random taxon pairs
    if abs(cor_pos) > abs(cor_non):

        boot_cors = []
        non_sig = False

        failed = 0
		
		# 1 error in 20 is the same as Pval = 0.05
        for b in range(20):
            temp = pd.concat([raw_numbers_pos.iloc[:, 0].sample(n=N, replace=True, axis=0).reset_index(drop=True),
                              raw_numbers_pos.iloc[:, 1].sample(n=N, replace=True, axis=0).reset_index(drop=True)],
                             axis=1)
            # boot_cors.append(temp.corr(method='pearson').fillna(0).iloc[1,0])
            if abs(temp.corr(method='pearson').fillna(0).iloc[1, 0]) > abs(cor_pos):
                failed += 1
                if failed >= 1:
                    non_sig = True
                    break
		
		# Record non-sign. coefficients just in case
        if not (non_sig):

            if abs(cor_pos) > abs(cor_non):
                print(gg)
                print(cor_pos)
                print(cor_non)
                with open('agp_hc_metab_%s.out' % sys.argv[1].split('_')[1][0], 'a') as f:
                    f.write("%s\t%s\t%s\t%s\t%s" % (gg[0], gg[1], cor_pos, cor_non, 'sig\n'))
        else:
            if abs(cor_pos) > abs(cor_non):
                with open('agp_hc_metab_%s.out' % sys.argv[1].split('_')[1][0], 'a') as f:
					f.write("%s\t%s\t%s\t%s\t%s" % (gg[0], gg[1], cor_pos, cor_non, 'non_sig\n'))
                print(gg)
                print('Non-signinficant')