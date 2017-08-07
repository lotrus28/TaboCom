GZ_PATH='/data5/bio/runs-galkin/PRJEB13679/hc'
ERR_PATH='/data5/bio/runs-galkin/PRJEB13679/ERRs'
FILTER_OUT='/data5/bio/runs-galkin/PRJEB13679/ERRs/filter_OUT'
BLOOM_PATH='/data5/bio/runs-galkin/PRJEB13679/qiime_16s/BLOOM.fasta'
PARAMS='/data5/bio/runs-galkin/PRJEB13679/qiime_16s/PARAMS'
GREENGENES_PATH='/data5/bio/runs-galkin/PRJEB13679/qiime_16s'
GREENGENES_OUT='/data5/bio/runs-galkin/PRJEB13679/greengenes_OUT'
OTU_TABLE_PATH='/data5/bio/runs-galkin/PRJEB13679/picked_otus'
COLLAPSED_TABLES_PATH='/data5/bio/runs-galkin/PRJEB13679/hc/collapsed'
PYTHON_SCRIPTS='/data5/bio/runs-galkin/PRJEB13679/qiime_16s'

# If data is preprocessed — no need to call taxons
# Just combine all the otu-tables and add pseudocounts
preprocessed=$1

# How to normalize data?
# 'div' for division
#'ln' for logarithmic
# 'none' for no normalisation
norm=$2

if [ $preprocessed = 'n' ]; then
	
	if [ ! -d $OTU_TABLE_PATH ]; then
		mkdir $OTU_TABLE_PATH
	fi
	
	if [ ! -d $ERR_PATH ]; then
		mkdir $ERR_PATH
	fi

	for FILE in $GZ_PATH/*.gz
	do
		# Unzip a downloaded fastq-file
		FILENAME=$(basename "${FILE}" .gz)
		if [ ! -f  $ERR_PATH/"${FILENAME}".fastq ]; then
			gunzip -c "${FILE}" > $ERR_PATH/"${FILENAME}"
			# Turn fastq to fasta
			# Disregard quality filtering
		fi
		FILENAME=$(basename "${FILE}" .fastq.gz)
		if [ ! -f $OTU_TABLE_PATH/$FILENAME ]; then
			awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $ERR_PATH/"${FILENAME}".fastq > $ERR_PATH/temp.fasta
			# Get rid of blooming bacteria
			echo "Picking OTUs for $FILENAME"
			sleep 60s
			pick_closed_reference_otus.py -i $ERR_PATH/temp.fasta -o $FILTER_OUT -r $BLOOM_PATH -p $PARAMS
			sleep 100s
			filter_fasta.py -f $ERR_PATH/temp.fasta -m $FILTER_OUT/uclust_ref_picked_otus/temp_otus.txt -n -o $ERR_PATH/bloom_filtered.fasta
			# Map filtered reads to greengenes DB
			pick_closed_reference_otus.py -i $ERR_PATH/bloom_filtered.fasta -o $GREENGENES_OUT -r $GREENGENES_PATH/gg_13_8_99.fasta -t $GREENGENES_PATH/gg_13_8_99.gg.tax -p $PARAMS -f
			FILENAME=$(basename "$FILE")
			echo "$FILENAME"
			sleep 90s
			mv $GREENGENES_OUT/otu_table.biom $OTU_TABLE_PATH/"${FILENAME%%.*}"
			echo "Made OTU table for $OTU_TABLE_PATH/${FILENAME%%.*}"
			
			until [ ! -d $GREENGENES_OUT ]
			do
				rm -rf $GREENGENES_OUT
			done
			sleep 5s
			rm -rf $FILTER_OUT
			sleep 5s
			rm $ERR_PATH/temp.fasta
			sleep 5s
			rm $ERR_PATH/bloom_filtered.fasta
			sleep 15s
		fi
	done
	
	if [ ! -d $COLLAPSED_TABLES_PATH ]; then
		mkdir $COLLAPSED_TABLES_PATH
	fi
	
	for FILE in $OTU_TABLE_PATH/*
	do
		FILENAME=$(basename "$FILE")
		if [ ! -f $COLLAPSED_TABLES_PATH/$FILENAME.txt ]; then
			python $PYTHON_SCRIPTS/collapse_biom.py $FILE $COLLAPSED_TABLES_PATH/$FILENAME.txt
		fi
	done
fi

if [ ! -d ./ADDED ]; then
	mkdir ./ADDED
fi

num_added=$(ls ./ADDED | wc -l)
num_to_add=$(ls $COLLAPSED_TABLES_PATH | wc -l)
if [  "$num_added" -ne "$num_to_add" ]; then
	for FILE in $COLLAPSED_TABLES_PATH/*
	do
		FILENAME=$(basename "$FILE" .txt)
		if [ ! -f ./ADDED/{$FILENAME}_added_higher.txt ]; then
			python $PYTHON_SCRIPTS/add_higher_taxons.py $FILE ./ADDED
		fi
	done
fi

# Some taxons may have several identifiers in greengenes
# So we need to combine all counts that have similar taxonomy but different IDs
if [ ! -f ./FINAL_OTU_TABLE.txt ]; then
	ls -d -1 ./ADDED/*.* > ./temp.txt
	echo 'Creating FINAL.txt'
	python $PYTHON_SCRIPTS/combine_collapsed_otu_tables.py ./temp.txt ./FINAL_OTU_TABLE.txt
 	rm ./temp.txt
fi

python $PYTHON_SCRIPTS/add_pseudocounts_normalize.py ./FINAL_OTU_TABLE.txt ./FINAL_transformed.txt $norm

