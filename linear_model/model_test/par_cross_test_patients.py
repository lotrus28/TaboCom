import pandas as pd
import itertools
import sys
import time
import multiprocessing
import numpy as np

num_compl = None

# params = 'E:\\Lab\\FHM\\AGP\\Parallel\\py\\PARAM_FINAL_OTU_TABLE_ibd.txt'
# ibd_folders = 'E:\\Lab\\FHM\\AGP\\Parallel\\FINAL_OTU_TABLE_ibd'
# heal_folders = 'E:\\Lab\\FHM\\AGP\\Parallel\\FINAL_OTU_TABLE_heal'
# out = 'E:\\Lab\\FHM\\AGP\\Parallel\\py\\OUT.txt'
# test_or_valid = 'test'

ibd_folders = sys.argv[1]
heal_folders = sys.argv[2]
out = sys.argv[3]
test_or_valid = sys.argv[4]

def init(args):
    global num_compl
    num_compl = args


def acc_for_multiplier(RMSEs, mult):
    heal_RMSE = RMSEs[0]
    ibd_RMSE = RMSEs[1]
    delta = heal_RMSE - (ibd_RMSE * mult)
    delta = ['Healthy' if i < 0 else 'IBD' for i in delta]
    return (delta)

def return_tax_names(ID_to_taxon, profiles):
    tax_code = {}
    with open(ID_to_taxon) as f:
        for line in f.read().splitlines():
            temp = line.split()
            # tax ID is in the 1st column tax name -- 2nd
            tax_code[temp[0]] = temp[1]
    profiles.index = [tax_code[x] for x in profiles.index]
    return (profiles)


def RMSE_calc(TEMP):
    i = TEMP[0]
    params = TEMP[1]

    global num_compl
    with num_compl.get_lock():
        num_compl.value += 1

    if (num_compl.value % 100 == 0):
        print( 'echo {} out of {}'.format(str(num_compl.value), str(params.shape[0])) )

    def profile_triple_SS(profile, models):
        SS = []
        for i in range(models.shape[0]):
            temp = models.iloc[i]
            try:
                pred = float(temp['Coef_p1']) * profile.loc[temp['Predictor1']] + \
                       float(temp['Coef_p2']) * profile.loc[temp['Predictor2']] + float(temp['Intercept'])
                SS.append((profile.loc[temp['Response']] - pred) ** 2)
            except KeyError:
                continue
        return (SS)

    def profile_pair_SS(profile, models):
        SS = []
        for i in models.index:
            temp = models.iloc[i]
            try:
                pred = temp['Coef_p'] * profile.loc[temp.loc['Predictor']] + temp['Intercept']
                SS.append((profile.loc[temp['Response']] - pred) ** 2)
            except KeyError:
                continue
        return (SS)

    def RMSE_test_profile(do_models_heal, tr_models_heal, do_models_ibd, tr_models_ibd, profile):
        heal_tr_ss = profile_triple_SS(profile, tr_models_heal)
        heal_do_ss = profile_pair_SS(profile, do_models_heal)
        heal_SS = list(heal_tr_ss) + list(heal_do_ss)
        heal_RMSE = (sum(heal_SS) / len(heal_SS)) ** (0.5)

        ibd_tr_ss = profile_triple_SS(profile, tr_models_ibd)
        ibd_do_ss = profile_pair_SS(profile, do_models_ibd)
        ibd_SS = list(ibd_tr_ss) + list(ibd_do_ss)
        ibd_RMSE = (sum(ibd_SS) / len(ibd_SS)) ** (0.5)

        # delta = heal_RMSE-(ibd_RMSE*multiplier)
        # delta = ['Healthy' if i <0 else 'IBD' for i in delta]
        return (heal_RMSE, ibd_RMSE)

    H_Tt = str(params.loc[i, 'H_Taxon_trim'])
    H_Sc = str(params.loc[i, 'H_SparCC_cor'])
    H_Sp = str(params.loc[i, 'H_SparCC_pval'])
    H_Ta = str(int(params.loc[i, 'H_Tax_adjacency']))
    H_Pf = str(params.loc[i, 'H_Pair_fstat'])
    #
    I_Tt = str(params.loc[i, 'I_Taxon_trim'])
    I_Sc = str(params.loc[i, 'I_SparCC_cor'])
    I_Sp = str(params.loc[i, 'I_SparCC_pval'])
    I_Ta = str(int(params.loc[i, 'I_Tax_adjacency']))
    I_Pf = str(params.loc[i, 'I_Pair_fstat'])

    # print('H_Tt: ' + str(H_Tt))
    # print('H_Sc: ' + str(H_Sc))
    # print('H_Sp: ' + str(H_Sp))
    # print('H_Ta: ' + str(H_Ta))
    # print('H_Pf: ' + str(H_Pf))

    # print('I_Tt: ' + str(I_Tt))
    # print('I_Sc: ' + str(I_Sc))
    # print('I_Sp: ' + str(I_Sp))
    # print('I_Ta: ' + str(I_Ta))
    # print('I_Pf: ' + str(I_Pf))

    template = '{}/trim_{}/pval_{}/dist_{}_cor_{}/fstat_{}'
    paths = {k: template.format(x, Tt, Sp, Ta, Sc, Pf)
             for k, x, Tt, Sp, Ta, Sc, Pf in
             (
                 ('ibd', ibd_folders, I_Tt, I_Sp, I_Ta, I_Sc, I_Pf),
                 ('heal', heal_folders, H_Tt, H_Sp, H_Ta, H_Sc, H_Pf)
             )}
    p_tax = {k: v for k, v in (
        ('ibd', ibd_folders + '/trim_%s/tax_code.txt' % I_Tt),
        ('heal', heal_folders + '/trim_%s/tax_code.txt' % H_Tt)
    )}
    p_do_models = {k: v for k, v in (
        ('ibd', paths['ibd'] + '/pair_models_{}.txt'.format(I_Pf)),
        ('heal', paths['heal'] + '/pair_models_{}.txt'.format(H_Pf))
    )}
    p_tr_models = {k: v for k, v in (
        ('ibd', paths['ibd'] + '/triplet_models_{}.txt'.format(I_Pf)),
        ('heal', paths['heal'] + '/triplet_models_{}.txt'.format(H_Pf))
    )}
    p_test_profile = {k: v for k, v in (
        ('ibd', ibd_folders + '/trim_{}/{}.txt'.format(I_Tt, test_or_valid)),
        ('heal', heal_folders + '/trim_{}/{}.txt'.format(H_Tt, test_or_valid))
    )}
    try:
        with open(p_tr_models['ibd']) as f:
            if f.readline() == 'No pair models to base triple models on':
                print('No pair models to base triple models on')
                params.loc[i, 'Sensitivity'] = -1
                params.loc[i, 'Specificity'] = -1
                sens = -1
                spec = -1
                return (i, sens, spec)
        with open(p_tr_models['heal']) as f:
            if f.readline() == 'No pair models to base triple models on':
                print('No pair models to base triple models on')
                params.loc[i, 'Sensitivity'] = -1
                params.loc[i, 'Specificity'] = -1
                sens = -1
                spec = -1
                return (i, sens, spec)

        do_models = {k: pd.read_table(v, header=0, index_col=None, sep='\t', engine="python") for k, v in
                     p_do_models.items()}
        tr_models = {k: pd.read_table(v, header=0, index_col=None, sep='\t', engine="python") for k, v in
                     p_tr_models.items()}
        test_profile = {k: pd.read_table(v, header=0, index_col=0, sep='\t', engine="python") for k, v in
                        p_test_profile.items()}

        for el in test_profile:
            test_profile[el] = return_tax_names(p_tax[el], test_profile[el])

        H = RMSE_test_profile(do_models['heal'], tr_models['heal'], do_models['ibd'], tr_models['ibd'],
                              test_profile['heal'])
        I = RMSE_test_profile(do_models['heal'], tr_models['heal'], do_models['ibd'], tr_models['ibd'],
                              test_profile['ibd'])


    except:
        sens = -1
        spec = -1
        return (i, sens, spec)

    multiplier = params.loc[i, 'RMSE_sign']

    H_RMSE = acc_for_multiplier(H, multiplier)
    I_RMSE = acc_for_multiplier(I, multiplier)

    sens = I_RMSE.count('IBD') / float(len(I_RMSE))
    spec = H_RMSE.count('Healthy') / float(len(H_RMSE))

    return (i, spec, sens)

# This script reads parameters form a txt file
# The file should contain a header for parameter combinations of both diseassed and helthy patients
try:
    params = pd.read_table('CROSS_PARAMETERS.txt', header=0, index_col=None, engine="python")
except:
    print("Calculating parameter table")
    columns = ['H_Taxon_trim', 'H_SparCC_pval', 'H_SparCC_cor', 'H_Tax_adjacency', 'H_Pair_fstat',
               'I_Taxon_trim', 'I_SparCC_pval', 'I_SparCC_cor', 'I_Tax_adjacency', 'I_Pair_fstat',
               'RMSE_sign', 'Specificity', 'Sensitivity']
    params = pd.DataFrame(columns=columns)

    values = {}

    values['H_Taxon_trim'] = np.r_[2:6] / 10.
    values['H_SparCC_pval'] = [0.001, 0.01, 0.05]
    values['H_SparCC_cor'] = np.r_[2:5] / 10.
    values['H_Tax_adjacency'] = np.r_[0:4]
    values['H_Pair_fstat'] = [0.001]
    #
    values['I_Taxon_trim'] = np.r_[2:6] / 10.
    values['I_SparCC_pval'] = [0.001, 0.01, 0.05]
    values['I_SparCC_cor'] = np.r_[2:5] / 10.
    values['I_Tax_adjacency'] = np.r_[0:4]
    values['I_Pair_fstat'] = [0.001]
    #
    values['RMSE_sign'] = [0.33, 0.5, 0.66, 1., 1.5, 2.0, 3.0]

    permutations = [x for x in itertools.product(
        values['H_Taxon_trim'],
        values['H_SparCC_pval'],
        values['H_SparCC_cor'],
        values['H_Tax_adjacency'],
        values['H_Pair_fstat'],
        values['I_Taxon_trim'],
        values['I_SparCC_pval'],
        values['I_SparCC_cor'],
        values['I_Tax_adjacency'],
        values['I_Pair_fstat'],
        values['RMSE_sign']
    )]

    # Make a file with all parameter combinations
    for set in permutations:
        newline = list(set) + [None, None]
        params = params.append(pd.Series(newline, index=columns), ignore_index=True)
    params.to_csv('CROSS_PARAMETERS.txt', sep='\t', index=False)

# Parallelize calculations with multiprocessing
if __name__ == '__main__':
    num_compl = multiprocessing.Value('i', 0)

    pool = multiprocessing.Pool(processes=25, initializer=init, initargs=(num_compl,))
    temp = ([row, params] for row in params.index)
    res = pool.map_async(RMSE_calc, temp, chunksize=1)
    res.wait()
    for el in res.get():
        i = el[0]
        spec = el[1]
        sens = el[2]
        params.loc[i, 'Sensitivity'] = sens
        params.loc[i, 'Specificity'] = spec
    params.to_csv(out, sep='\t', index=False)
