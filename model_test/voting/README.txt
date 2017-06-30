= voting_test.R
A script for voting based model validation.
An error is calculated for each model and sample. This error is compared to errors, produced while applying models to training data.
P-values are calculated based on how many errors in training data are bigger than the error of the sample. eg: if a helthy cohort-derived model has an extremely small P-value for a sample (error is significantly bigger than expected for this model), then the sample is less likely to be from healthy cohort.
Each model casts a vote based on P-value threshold. eg: if >50% models derived from healthy cohort data have small errors — the sample is healthy.

= CROSS_PARAMETERS.txt
List of parameters assessed in FINAL_OTU_TABLE_heal and FINAL_OTU_TABLE_ibd

=FINAL_OTU_TABLE_*
The subfolders contain training and testing datasets as well as models