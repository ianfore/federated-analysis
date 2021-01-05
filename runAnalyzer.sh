#!/bin/bash

if [ $# -ne 1 ]
then
	echo "usage $0 chr-num"
	exit 1
fi
CHR=$1
#NUM_CORES=$(grep -c processor /proc/cpuinfo)
NUM_CORES=1

DATA_DIR=/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data

python app/cooccurrence/cooccurrenceAnalyzer.py \
-o $DATA_DIR/$CHR/$CHR-out.json \
-v $DATA_DIR/$CHR/$CHR-vpi.json \
-i $DATA_DIR/$CHR/$CHR-ipv.json \
-b ./data/brca-variants.tsv \
-r ./data/brca-regions.json \
-a $DATA_DIR/ancestries.json \
-d $DATA_DIR/$CHR/OUTPUT \
-n $NUM_CORES \
