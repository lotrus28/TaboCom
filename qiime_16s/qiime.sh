GZ_PATH='/data5/bio/runs-galkin/AGP/QIIME/test1'
ERR_PATH='/data5/bio/runs-galkin/AGP/QIIME/ERRs'
FILTER_OUT='/data5/bio/runs-galkin/AGP/QIIME/ERRs/filter_OUT'
BLOOM_PATH='/data5/bio/runs-galkin/AGP/QIIME/BLOOM.fasta'
PARAMS='/data5/bio/runs-galkin/AGP/QIIME/PARAMS'
GREENGENES_PATH='/data5/bio/runs-galkin/AGP/QIIME'
GREENGENES_OUT='/data5/bio/runs-galkin/AGP/QIIME/greengenes_OUT'
OTU_TABLE_PATH='/data5/bio/runs-galkin/AGP/QIIME/picked_otus'
COLLAPSED_TABLES_PATH='/data5/bio/runs-galkin/AGP/QIIME/heal/collapsed'
PYTHON_SCRIPTS='/data5/bio/runs-galkin/AGP/QIIME/py'

# If data is preprocessed — no need to call taxons
# Just combine all the otu-tables and add pseudocounts
preprocessed=$1

# How to normalize data?
# 'div' for division
#'ln' for logarithmic
# 'none' for no normalisation
norm=$2

if [ preprocessed = 'n' ]; then
	mkdir $OTU_TABLE_PATH
	mkdir $ERR_PATH

	for FILE in $GZ_PATH/*.gz
	do
		# Unzip a downloaded fastq-file
		FILENAME=$(basename "${FILE}" .gz)
		gunzip -c "${FILE}" > $ERR_PATH/"${FILENAME}"
		# Turn fastq to fasta
		# Disregard quality filtering
		awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $ERR_PATH/"${FILENAME}" > $ERR_PATH/temp.fasta
		# Get rid of blooming bacteria
		pick_closed_reference_otus.py -i $ERR_PATH/temp.fasta -o $FILTER_OUT -r $BLOOM_PATH -p $PARAMS
		filter_fasta.py -f $ERR_PATH/temp.fasta -m $FILTER_OUT/uclust_ref_picked_otus/temp_otus.txt -n -o $ERR_PATH/bloom_filtered.fasta
		rm -rf $FILTER_OUT
		rm $ERR_PATH/temp.fasta
		# Map filtered reads to greengenes DB
		pick_closed_reference_otus.py -i $ERR_PATH/bloom_filtered.fasta -o $GREENGENES_OUT -r $GREENGENES_PATH/gg_13_8_99.fasta -t $GREENGENES_PATH/gg_13_8_99.gg.tax -p $PARAMS
		FILENAME=$(basename "$FILE")
		echo "$FILENAME"
		mv $GREENGENES_OUT/otu_table.biom $OTU_TABLE_PATH/"${FILENAME%%.*}"
		echo "$OTU_TABLE_PATH/${FILENAME%%.*}"
		rm -rf $GREENGENES_OUT
		rm $ERR_PATH/bloom_filtered.fasta
	done

	mkdir $COLLAPSED_TABLES_PATH
	for FILE in $OTU_TABLE_PATH/*
	do
		FILENAME=$(basename "$FILE")
		python3 $PYTHON_SCRIPTS/collapse_biom.py $FILE $COLLAPSED_TABLES_PATH/$FILENAME.txt
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
		python3 $PYTHON_SCRIPTS/add_higher_taxons.py $FILE ./ADDED
	done
fi

# Some taxons may have several identifiers in greengenes
# So we need to combine all counts that have similar taxonomy but different IDs
if [ ! -f ./FINAL_OTU_TABLE.txt ]; then
	ls -d -1 ./ADDED/*.* > ./temp.txt
	python3 $PYTHON_SCRIPTS/combine_collapsed_otu_tables.py ./temp.txt ./FINAL_OTU_TABLE.txt
	rm ./temp.txt
fi

python3 $PYTHON_SCRIPTS/add_pseudocounts_normalize.py ./FINAL_OTU_TABLE.txt ./FINAL_transformed.txt $norm

