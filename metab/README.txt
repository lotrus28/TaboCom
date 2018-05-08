Thåse scripts perform metabolic complementarity analysis on correlation table for taxa.

== ./reference
All the metabolic reference files are derived from nature.com/articles/ncomms15393.
The folder includes:
	- group_to_species.txt: a dict of metabolic group names (Acetoge, Methanogen etc) to species present in the reference;
	- compound_to_species.txt: a table of compound producing/consuming capabilities for reference species;
	- prod.txt: a table tof size NxM, where N is the # of taxa and M is the # of compounds. Each cell contains the probability that a random species from the corresponding taxa can produce a metabolite;
	- consum.txt: -//- but for consumer species

== sequential_metabolic_corrs.py
The script to sequentially calculate metabolic complementarity (MC) among all compounds.
In: consument and producent reference tables; taxon correlation and a corresponding significance tables;
Out: *.out file containing metabolic complementarity coef-s that are deemed significant @ 0.05 level

== parallel_metabo_corrs_master.py
The script that launches MC calculation parallelised between several nodes on sun grid engine
-//-inputs
Creates files containing compounds each node has to calculate MC for.

== parallel_metabo_corrs_subordinate.py
The script that calculates MC for compounds specified in steamN.txt files.
Outputs are a number of tables containing MC for different compounds.

== pval_metab.py
The script that computes exact Pvalue for MC coeffficients
Worsk the same way as `parallel_metabo_corrs_subordinate.py`.