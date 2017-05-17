= process.py
This script takes in taxonomy tree used at KEGG (same as NCBI Taxonomy) and predicitive models files for two cohorts.
The output is 'taxon_kegg.csv' that contains info on the fullness of metabolic pathways in taxons, present in models.
This script is not that useful

= taxon_kegg.csv
Example output of 'process.py'

= HEALTHY_pair_models_0.001.txt
= HEALTHY_triplet_models_0.001.txt
= IBD_pair_models_0.001.txt
= IBD_triplet_models_0.001.txt
Example input to 'process.py'

= kegg.py
This script contains functions that calculate metabolic vectors for given pathways and taxons.
It also can assess metabolic distance difference between model-associated taxons and non-model taxons. Significant differences are printed output

= *_models.txt
= *_tax_code.txt
= *_train.txt
Sample inputs to 'kegg.py'

= br08610.keg
= kos.txt
= gene_links.txt
Data parsed from KEGG.
The scripts are modified to pull data from KEGG if these files are missing

= ko00780.txt
= IBD_ko00780.txt
Sample outputs containing normalized metabolic vectors for IBD-model associated taxons.

Disclaimer:
I fully understand these scripts are undercommrnted and need more clarity.
We're working at it
