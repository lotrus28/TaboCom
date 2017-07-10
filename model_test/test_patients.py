import pandas as pd
import numpy as np
import sys

def profile_triple_SS(profile, models):
    SS = []
    for i in range(models.shape[0]):
        temp = models.iloc[i]
        try:
            pred = float(temp['Coef_p1'])*profile.loc[temp['Predictor1']] + \
                   float(temp['Coef_p2'])*profile.loc[temp['Predictor2']] + float(temp['Intercept'])
            SS.append((profile.loc[temp['Response']] - pred) ** 2)
        except KeyError:
            continue
    return(SS)

def profile_pair_SS(profile, models):
    SS = []
    for i in models.index:
        temp = models.iloc[i]
        try:
            pred = temp['Coef_p'] * profile.loc[temp.loc['Predictor']] + temp['Intercept']
            SS.append((profile.loc[temp['Response']] - pred) ** 2)
        except KeyError:
            continue
    return(SS)

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
    return(heal_RMSE, ibd_RMSE)

def acc_for_multiplier(RMSEs, mult):
    heal_RMSE = RMSEs[0]
    ibd_RMSE = RMSEs[1]
    delta = heal_RMSE - (ibd_RMSE * mult)
    delta = ['Healthy' if i < 0 else 'IBD' for i in delta]
    return(delta)

def return_tax_names(ID_to_taxon, profiles):
    tax_code = {}
    with open(ID_to_taxon) as f:
        for line in f.read().splitlines():
            temp = line.split()
            # tax ID is in the 1st column tax name -- 2nd
            tax_code[temp[0]] = temp[1]
    profiles.index = [tax_code[x] for x in profiles.index]
    return(profiles)

# params = 'E:\\Lab\\FHM\\AGP\\Parallel\\py\\PARAM_FINAL_OTU_TABLE_ibd.txt'
# ibd_folders = 'E:\\Lab\\FHM\\AGP\\Parallel\\FINAL_OTU_TABLE_ibd'
# heal_folders = 'E:\\Lab\\FHM\\AGP\\Parallel\\FINAL_OTU_TABLE_heal'

params = sys.argv[1]
ibd_folders = sys.argv[2]
heal_folders = sys.argv[3]
out = sys.argv[4]
test_or_valid = sys.argv[5]

params = pd.read_table(params, header = 0, index_col = None)
print(len(params.index))
Pf_prev = -1
for i in params.index:

    Tt = str(params.loc[i,'Taxon_trim'])
    Sc = str(params.loc[i,'SparCC_cor'])
    Sp = str(params.loc[i,'SparCC_pval'])
    Ta = str(int(params.loc[i,'Tax_adjacency']))
    Pf = str(params.loc[i,'Pair_fstat'])

    if Pf_prev != Pf:

        Pf_prev = Pf

        template = '{}/trim_{}/pval_{}/dist_{}_cor_{}/fstat_{}'
        paths = {k:template.format(x, Tt, Sp, Ta, Sc, Pf)
                 for k,x in (('ibd',ibd_folders), ('heal',heal_folders))}
        print(paths)
        p_tax = {k:v for k,v in (
                                ('ibd',ibd_folders+'/trim_%s/tax_code.txt' % Tt),
                                ('heal',heal_folders+'/trim_%s/tax_code.txt' % Tt)
                                )}
        p_do_models = {k:v for k,v in (
                                        ('ibd', paths['ibd'] +'/pair_models_{}.txt'.format(Pf)),
                                        ('heal', paths['heal'] + '/pair_models_{}.txt'.format(Pf))
                                     )}
        p_tr_models = {k:v for k,v in (
                                        ('ibd', paths['ibd'] + '/triplet_models_{}.txt'.format(Pf)),
                                        ('heal',paths['heal'] + '/triplet_models_{}.txt'.format(Pf))
                                    )}
        p_test_profile = {k:v for k,v in (
                                        ('ibd', ibd_folders + '/trim_{}/{}.txt'.format(Tt,test_or_valid)),
                                        ('heal', heal_folders + '/trim_{}/{}.txt'.format(Tt,test_or_valid))
                                        )}
        try:
            with open(p_tr_models['ibd']) as f:
                if f.readline() == 'No pair models to base triple models on':
                    print('No pair models to base triple models on')
                    params.loc[i, 'Sensitivity'] = 'NA'
                    params.loc[i, 'Specificity'] = 'NA'
                    continue
            with open(p_tr_models['heal']) as f:
                if f.readline() == 'No pair models to base triple models on':
                    print('No pair models to base triple models on')
                    params.loc[i, 'Sensitivity'] = 'NA'
                    params.loc[i, 'Specificity'] = 'NA'
                    continue

            do_models = {k:pd.read_table(v, header=0, index_col=None, sep='\t') for k,v in p_do_models.items()}
            tr_models = {k:pd.read_table(v, header=0, index_col=None, sep='\t') for k, v in p_tr_models.items()}
            test_profile = {k:pd.read_table(v, header=0, index_col=0, sep='\t') for k,v in p_test_profile.items()}
            # tax = {k:pd.read_table(v, header=None, index_col=0, sep='\t') for k,v in tax.items()}

            for el in test_profile:
                test_profile[el] = return_tax_names(p_tax[el], test_profile[el])

            H = RMSE_test_profile(do_models['heal'], tr_models['heal'], do_models['ibd'], tr_models['ibd'],
                                  test_profile['heal'])
            I = RMSE_test_profile(do_models['heal'], tr_models['heal'], do_models['ibd'], tr_models['ibd'],
                                  test_profile['ibd'])

        except:
            print('File misssing')
            params.loc[i, 'Sensitivity'] = 'Miss'
            params.loc[i, 'Specificity'] = 'Miss'
            continue

    for file in (p_tr_models['ibd'], p_tr_models['heal']):
        try:
            with open(file) as f:
                if f.readline() == 'No pair models to base triple models on':
                    print('No pair models to base triple models on')
                    params.loc[i, 'Sensitivity'] = 'NA'
                    params.loc[i, 'Specificity'] = 'NA'
                    continue
        except:
            print('File misssing')
            params.loc[i, 'Sensitivity'] = 'Miss'
            params.loc[i, 'Specificity'] = 'Miss'
            continue


    multiplier = params.loc[i,'RMSE_sign']

    H_RMSE = acc_for_multiplier(H, multiplier)
    I_RMSE = acc_for_multiplier(I, multiplier)

    sens = I_RMSE.count('IBD') / float(len(I_RMSE))
    spec = H_RMSE.count('Healthy') / float(len(H_RMSE))
    params.loc[i, 'Sensitivity'] = sens
    params.loc[i, 'Specificity'] = spec
    print( str(i) + ' out of ' + str(params.shape[0]) )

params.to_csv(out, sep='\t', index=False)