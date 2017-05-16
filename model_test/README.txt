This folder contains scripts required for testing predictive model accuracy

= test_patients.py
This scripts takes in a parameter table and tests predictive models obtained with corresponding parameters.
The script calculates specificity and sensitivity of test.
Predictive models for both cohorts are bound to be produced through the same parameters.
You can specify whether to use validation set or testing set with the last parameter set to 'valid' or 'test'
The script ignores missing model files
Script launch command:
python3 test_patients.py ./PARAM_FINAL_OTU_TABLE_ibd.txt FINAL_OTU_TABLE_ibd FINAL_OTU_TABLE_heal test.txt test

= par_cross_test_patients.py
The script requires a 'CROSS_PARAMETERS.txt' to work. 'CROSS_PARAMETERS.txt' is a table of parameter combinations to use while testing. Unlike 'test_patients.py' this script uses predictive models from two cohorts that were obtained with different parameters. E.g. it can use models with P-value = 0.01 for healthy cohort and P-value = 0.001 for IBD-cohort, while 'test_patients.py' would use the same P-value for both cohorts.
The script produces this file if it is missing. In this case pay extra attention to parameter combinations that will be present in the resulting table.
Script launch command:
python3 par_cross_test_patients.py /data5/bio/runs-galkin/Parallel/FINAL_OTU_TABLE_ibd /data5/bio/runs-galkin/Parallel/FINAL_OTU_TABLE_heal ./test.txt test

= test.txt
This is a sample output of test_patients.py

= PARAM_FINAL_OTU_TABLE_ibd.txt
This is a sample parameter combinations set