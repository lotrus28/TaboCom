= make_err_table.R
A script to calculate prediction errors of previosly picked linear models.

= boruta_test.R
A script to discover descriptive models (as relevant fetures in rabdom forest estimation) and create a classifier based on these models.

= ./input_data/
All tables needed to build a classifier (otu tables, all possible models calculated etc)

= ./forest_related/
Outputs of make_err_table.R ('errs_to_boruta.txt') and boruta_test.R ('randomForest.Rdata')
