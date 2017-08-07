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
values['Taxon_trim'] =[0.5]
values['SparCC_cor'] = [0.2]
values['SparCC_pval'] = [0.05]

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
    teach_predictor(path, [i, min(values['SparCC_pval']), 100], [0,1,1], job, None)
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
    teach_predictor(path, [i[0], i[1], 0.2], [1,0,1], job, wait)
    time.sleep(5)
wait = job

# time.sleep(3600)

# And we go all the way down
job = "Filter_sig"
Tt_Sp_Sc_Ta = [x for x in itertools.product(values['Taxon_trim'],
                                            values['SparCC_pval'],
                                            values['SparCC_cor'])]
for i in Tt_Sp_Sc_Ta:
    teach_predictor(path, [i[0], i[1], i[2]], [1,1,0], job, wait)
    time.sleep(3)
wait = job

time.sleep(1000)

wait = job
for f in ["Rplots.pdf", "cov_mat_SparCC.out", "get_pair_fstats.Rout"]:
    command = 'echo "rm %s" | qsub -N CleanUp -hold_jid pair_Fstat -cwd' % f
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
command = 'echo "find . -maxdepth 1 -type f -size 0 | xargs -d"\n" rm -f" | qsub -N Finish -hold_jid CleanUp -cwd'
subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

# ls | grep -P '\.o' | xargs -d"\n" rm -f
# find . -maxdepth 1  -type f -size 0 | xargs -d"\n" rm -f
