import itertools
import os.path
import sys
import subprocess
import time
import fileinput

import numpy as np
import pandas as pd

# Enter 1 parameter: otu table with reads
path = sys.argv[1]
cond = path.split('/')[-1].split('.')[0]

def teach_predictor(path, params, same, job, wait):
    time.sleep(1)
    if not(wait is None):
        wait = ' -hold_jid ' + wait + ' -cwd'
        command = 'echo "bash ./teach_models.sh ' + path + ' ' + ' '.join([str(x) for x in params]) + ' ' + \
                  ' '.join([str(x) for x in same]) +'" | qsub -N ' + job + wait
    else:
        command = 'echo "bash ./teach_models.sh ' + path + ' ' + ' '.join([str(x) for x in params]) + ' ' + \
                  ' '.join([str(x) for x in same]) +  '" | qsub -N ' + job + ' -cwd'
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # input_string = "./teach_models.sh {} {} {}\n".format(*(map(pipes.quote, [path, params, same])))
    # print(input_string)
    # p = subprocess.Popen(['qsub', '-N', job, '-cwd'], stdin=subprocess.PIPE)
    # out, err = p.communicate(input_string)
    return()

# Create a df with all possible combinations of parameters
columns = ['Taxon_trim', 'SparCC_pval', 'SparCC_cor', 'Tax_adjacency',
           'Pair_fstat', 'RMSE_sign', 'Specificity', 'Sensitivity']
params = pd.DataFrame(columns=columns)

values = {}
values['Taxon_trim'] = np.r_[2:6] / 10.
values['SparCC_cor'] = np.r_[2:5] / 10.
values['SparCC_pval'] = [0.001, 0.01, 0.05]
values['Tax_adjacency'] = np.r_[0:4]
values['Pair_fstat'] = [0.001]
values['RMSE_sign'] = [0.33, 0.5, 0.66, 1., 1.5, 2.0, 3.0]

# values['Taxon_trim'] = [0.3]
# values['SparCC_cor'] = [0.3]
# values['SparCC_pval'] = [0.001,0.01]
# values['Tax_adjacency'] = [2,3]
# values['Pair_fstat'] = [0.01,0.001]
# values['RMSE_sign'] = [1.]


# Make a parameter table in cases there is none
# It will be needed for test_patients.py
filename = 'PARAM_' + cond + '.txt'
if not(os.path.exists(filename)):
    permutations = [x for x in itertools.product(values['Taxon_trim'],
                                                values['SparCC_pval'],
                                                values['SparCC_cor'],
                                                values['Tax_adjacency'],
                                                values['Pair_fstat'],
                                                 values['RMSE_sign'])
                    ]
    # Make a file with all parameter combinations
    for set in permutations:
        newline = list(set) + [None, None]
        params = params.append(pd.Series(newline, index=columns), ignore_index=True)
    params.to_csv(filename, sep='\t', index=False)

# Erase brackets from tax_code, because later model calling can't work with them
# sed might not work for some reason

command = "sed -ie 's/\[/BRA/g;s/\]/KET/g;s/-/SLASH/g' %s" % path
subprocess.Popen(command, shell=True)
time.sleep(5)
command = "rm %se" % path
subprocess.Popen(command, shell=True)

with fileinput.FileInput(path, inplace=True, backup='.bak') as file:
    for line in file:
        print(line.replace(']', 'KET').replace('[','BRA').replace('-','SLASH'), end='')

# First, create SparCC-outputs
job = "Tax_trim"
for i in values['Taxon_trim']:
    teach_predictor(path, [i, min(values['SparCC_pval']), 100, 100, 100], [0,1,1,1,1], job, None)
wait = job

# time.sleep(10)

# Then one level lower
job = "Spar_pval"
Tt_Sp = [x for x in itertools.product(values['Taxon_trim'],
                                        values['SparCC_pval'])]
# Make it so first all SparCC with lowest P-value are calculated
# Otherwise
Tt_Sp = sorted(Tt_Sp, key=lambda x: x[1])
for i in Tt_Sp:
    teach_predictor(path, [i[0], i[1], 100, 100, 100], [1,0,1,1,1], job, wait)
    time.sleep(5)
wait = job

time.sleep(3600)

# And we go all the way down
job = "Filter_sig"
Tt_Sp_Sc_Ta = [x for x in itertools.product(values['Taxon_trim'],
                                            values['SparCC_pval'],
                                            values['SparCC_cor'],
                                            values['Tax_adjacency'])]
for i in Tt_Sp_Sc_Ta:
    teach_predictor(path, [i[0], i[1], i[2], int(i[3]), -1], [1,1,0,0,1], job, wait)
    time.sleep(3)
wait = job

time.sleep(1000)

job = 'Calc_models'
for i in values['Taxon_trim']:
    fstat_min = str(min(values['SparCC_pval']))
    fstat_max = str(max(values['SparCC_pval']))
    cor_thr = str(min(values['SparCC_cor']))
    p_sigcor = './{}/trim_{}/pval_{}/sig_cor_{}.txt'.format(cond,str(i),fstat_max,fstat_max)
    p_counts = './{}/trim_{}/train.txt'.format(cond,str(i))
    p_taxcode = './{}/trim_{}/tax_code.txt'.format(cond,str(i))
    out = './{}/trim_{}/'.format(cond,str(i))

    command = '''echo 'R CMD BATCH "--args {} {} {} {} {} {} {}" ./Scripts/calculate_all.models.R' '''.format(p_sigcor, p_counts, p_taxcode, fstat_min, fstat_max, cor_thr, out)
    command = command + '| qsub -N ' + job + ' -hold_jid ' + wait + ' -cwd'
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    time.sleep(5)

time.sleep(3600)
wait = job

job = 'pair_Fstat'
Tt_Sp_Sc_Ta_Pf = [x for x in itertools.product(values['Taxon_trim'],
                                                values['SparCC_pval'],
                                                values['SparCC_cor'],
                                                values['Tax_adjacency'],
                                                values['Pair_fstat'])]
for i in Tt_Sp_Sc_Ta_Pf:
    teach_predictor(path, [i[0], i[1], i[2], i[3], i[4]], [1,1,1,1,0], job, wait)
    time.sleep(1)

for f in ["Rplots.pdf", "cov_mat_SparCC.out", "get_pair_fstats.Rout"]:
    command = 'echo "rm %s" | qsub -N CleanUp -hold_jid pair_Fstat -cwd' % f
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
command = 'echo "find . -maxdepth 1 -type f -size 0 | xargs -d"\n" rm -f" | qsub -N Finish -hold_jid CleanUp -cwd'
subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

# ls | grep -P '\.o' | xargs -d"\n" rm -f
# find . -maxdepth 1  -type f -size 0 | xargs -d"\n" rm -f
