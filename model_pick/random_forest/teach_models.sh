#!/bin/bash
SCRIPT_PATH='/data5/bio/runs-galkin/Parallel/randScripts'
SPARCC_PATH='/data5/bio/runs-galkin/Parallel/SparCC_source'

collapsed_otu_table=$1

# Parameter values
# They are used to create and adress folders
Taxon_trim=$2
Sparcc_pval=$3
Sparcc_cor=$4

input_folder="$(dirname $collapsed_otu_table)"
output_subfolder=${collapsed_otu_table##*/}
output_subfolder=${output_subfolder%.*}
out_all="$input_folder/$output_subfolder"
out_trim="$out_all/trim_$Taxon_trim"
out_SparCC="$out_trim/Results_SparCC"


# Parameters being the same
# 1 or 0
Tt=$5
Sp=$6
Sc=$7

ITER=10
BOOT_ITER=$(echo 1/$Sparcc_pval | bc)
last=$(expr $BOOT_ITER - 1)
job_prefix=$Taxon_trim

echo "$Taxon_trim $Sparcc_pval $Sparcc_cor $Tax_adjacency $Pair_fstat"
echo "$Tt $Sp $Sc $Ta $Pf"

# If no precalculations have been made
# then make a new output folder
if [ "$Tt" -eq 0 ]; then

	if [ ! -d $out_all/ ]; then
		mkdir $out_all/
	fi
	
	if [ ! -d $out_trim/ ]; then
		mkdir $out_trim/
		# Produces 3 otu tables: 1 w.training set patient, 2nd w.test set and 3d w.validation set
		# For convenience taxonomies are reduced to short tags. Tag-taxonomy correspondence is written to'tax_code.txt'
		python3 $SCRIPT_PATH/separate_sets_and_shorten_tax_names.py $collapsed_otu_table $out_trim $Taxon_trim
		# R CMD BATCH "--args train.txt  $out_trim" $SCRIPT_PATH/deseq_norm.R
		# echo "Data cut and normalized"
		sleep 20s
	fi
	
	if [ ! -d $out_SparCC ]; then
		mkdir $out_SparCC/
		mkdir $out_SparCC/Resamplings/
		mkdir $out_SparCC/Bootstraps/
	fi
	
	if [ ! -f $out_SparCC/cor_mat_SparCC.out ]; then
		# SparCC is used to get correlation table and pvalues of correlation significance
		inp_SparCC="$out_trim/train.txt"
		if [ ! -f $out_SparCC/cor_mat_SparCC.out ]; then
			python $SPARCC_PATH/SparCC.py $inp_SparCC -i $ITER --cor_file=$out_SparCC/cor_mat_SparCC.out > $out_SparCC/sparcc.log
		fi
		
		if [ ! -f $out_SparCC/Resamplings/resamp_$last.txt ]; then
			python $SPARCC_PATH/MakeBootstraps.py $inp_SparCC -n $BOOT_ITER -t resamp_#.txt -p $out_SparCC/Resamplings/ >> $out_SparCC/sparcc.log
		fi
	fi
	
	last=$(($BOOT_ITER - 1))
	num_straps=$(ls ${out_SparCC}/Bootstraps | wc -l)
	if [ $num_straps -lt $BOOT_ITER ]; then
	
		declare -A hostnames
		hostnames+=( [0]="node6.net0.pyxis.ripcm.com" [1]="node8.net0.pyxis.ripcm.com" [2]="node9.net0.pyxis.ripcm.com" [3]="node10.net0.pyxis.ripcm.com" [4]="node11.net0.pyxis.ripcm.com" [5]="node12.net0.pyxis.ripcm.com" )
		for ((j = 0; j <= $(($BOOT_ITER/6)); j++))
		do
			for ((i = 0; i < 6; i++))
			do
				taskno=$(($j*6+$i))
				if [ $taskno -eq $BOOT_ITER ]; then
					break
				fi
				if [ ! -f $out_SparCC/Bootstraps/sim_cor_$taskno.txt ]; then
					tasknoW=$(($j*6+$i-6))
					node=${hostnames[$i]}
					echo "python $SPARCC_PATH/SparCC.py $out_SparCC/Resamplings/resamp_$taskno.txt -i $ITER --cor_file=$out_SparCC/Bootstraps/sim_cor_$taskno.txt >> $out_SparCC/sparcc.log" | qsub -pe make 3 -N job_${job_prefix}_${taskno} -hold_jid job_${job_prefix}_${tasknoW} -cwd -l hostname=$node -r yes
				fi
			done
		done
	fi
	
	boot_samples=$(ls -1 $out_SparCC/Bootstraps | wc -l)
	if [ $boot_samples -lt $BOOT_ITER ]; then
		declare -A reserved_hostnames
		reserved_hostnames+=( [0]="node1.net0.pyxis.ripcm.com" [1]="node2.net0.pyxis.ripcm.com" [2]="node4.net0.pyxis.ripcm.com" [3]="node5.net0.pyxis.ripcm.com")
		until [ $boot_samples -ge $BOOT_ITER ]
		do
			echo "$boot_samples bootstraps are ready out of $BOOT_ITER"
			sleep 600s
			boot_samples_new=$(ls -1 $out_SparCC/Bootstraps | wc -l)
			if [ $boot_samples_new -eq $boot_samples ]; then
				N_pending_jobs=$(qstat | grep '\sqw' | grep 'job_' | wc -l)
				pending_jobs=$(qstat | grep '\sqw' | grep 'job_' | awk '{print $3}')
				if [ $N_pending_jobs -eq 0 ];then
					N_pending_jobs=$(qstat | grep 'qw' | grep 'job_' | wc -l)
					pending_jobs=$(qstat | grep 'qw' | grep 'job_' | awk '{print $3}')
				fi
				if [ ! $N_pending_jobs -eq 0 ]; then
						for ((k = 1; k <= $N_pending_jobs; k++))
						do
							pending_job=$(echo $pending_jobs | awk '{print $k}')
							node=${reserved_hostnames[$((k%4))]}
							qalter -l hostname=$node -hold_jid no_wait $pending_job
							sleep 1s
						done
						echo "Resubmitted $N_pending_jobs jobs"
					fi
			fi
			boot_samples=$boot_samples_new
		done

	fi
	
	if [ ! -f $out_SparCC/pvals_two_sided.txt ] && [ -d $out_SparCC ]; then
		
		echo "python $SPARCC_PATH/PseudoPvals.py $out_SparCC/cor_mat_SparCC.out $out_SparCC/Bootstraps/sim_cor_#.txt $BOOT_ITER -o $out_SparCC/pvals_two_sided.txt -t 'two_sided'  >> $out_SparCC/sparcc.log" | qsub -N Pval_calc_${Taxon_trim}_${Sparcc_pval} -cwd
		
		# if [ -f $out_SparCC/sparcc.log ]; then
			# rm -rf $out_SparCC/Resamplings
			# rm -rf $out_SparCC/Bootstraps
			# rm $out_SparCC/sparcc.log
		# fi
		
	fi
fi
		
if [ "$Sp" -eq 0 ]; then
	# I cant think of a better way to wait for pval calculations finish
	# Somehow hold_jid fails
	until [ -f  $out_SparCC/pvals_two_sided.txt ]
	do
		sleep 30s
	done
	sleep 30s
	
	# We need to filter out insignificant correlations w.R-script
	if [ ! -f $out_trim/sig_cor_$Sparcc_pval.txt ]; then
		R CMD BATCH "--args $out_SparCC/pvals_two_sided.txt $out_SparCC/cor_mat_SparCC.out $out_trim/tax_code.txt $out_trim $Sparcc_pval" $SCRIPT_PATH/sig_cor.R
	fi

	if [ ! -f $out_trim/all_models.txt ]; then
		R CMD BATCH "--args $out_trim/sig_cor_$Sparcc_pval.txt $out_trim/train.txt $out_trim/tax_code.txt $Sparcc_cor $out_trim/" $SCRIPT_PATH/calculate_all.models.R
	fi
	
	# Let's see if there any more descriptive triple linear models
	if [ -f $out_trim/all_models.txt ]; then
		if [ ! -f $out_trim/all_triplet_models.txt ]; then
			python3 $SCRIPT_PATH/triplets_wei.py $out_trim/all_models.txt $out_trim/train.txt $out_trim/tax_code.txt $out_trim
		fi
	fi
	
fi
exit
