This folder contains scripts that generate predictive models.

= teach_models.sh
This script calculates predictive models for a set number of parameters
Script launch command:
./teach_models.sh FINAL_transformed_heal_ln.txt 0.3 0.01 0.3 3 0.001 0 0 0 0 0
Sun Grid Engine is required to run teach_models.sh as it carries out some calculations in parallel.
First 5 numbers stand for parameters to obtain predictive model.
Last 5 numbers stand for calculations required. In case you have the output that utilizes the 1st parameter you can replace '0 0 0 0 0' with '1 0 0 0 0' in orderr to avoid repeat calculations.
Pay attention to hostnames in this script: you can modify them to obtain the optimal resource usage.
The script is modified to use DeSeq-normalized counts. You can use raw count by commenting the next line in script:
R CMD BATCH "--args train.txt  $out_trim" $SCRIPT_PATH/deseq_norm.R

= FINAL_OTU_TABLE_heal.txt
= FINAL_OTU_TABLE_ibd.txt
These are sample inputs.
They are generated from sample fastq files via qiime_16s workflow

= parallel_script.py
This is a master script, that provides you the opportunity to calculate predictive models for multiple parameter combinations. 
It also creates a parameter table that is required for further model testing.
The output is a tree-like directory structure with predictive models in leaf folders
The script is configured to check for existing output. In case of error or code arrest you can restart it and the script will continue working after checking for existing output.
Pay attention to parameter combinations hardcoded in this script. You can modify them to calculate fewer or more predictive models.
Script launch command:
python3 parallel_script.py FINAL_OTU_TABLE_heal.txt

= Scripts
These contain various scripts that are launched by teach_models.sh and that carry out predictive model calculations.

= SparCC_source
This folder contains SparCC tool developed elsewhere that permits compositional data correlation analysis.
Several scripts have been modified to avoid errors specific to our laboratory server. In case these scripts fail to work feel free to utilize the original ones.

=scrambled_samples
This folder contains a script to generate scrambled testing samples to validate our model.
This script considers total sample to have more than 40 samples as it does not scramble the last 40 samples in each OTU-table.
These 40 samples will later be set aside as a testing set.