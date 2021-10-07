#!/bin/bash


DATA_DIR=/private/groups/patenlab/jcasalet/TOPMED_WORKDIR/DATA/F8/GT40/13
CHROM=13
GENE=BRCA2
P2=0.001
VCF_FILENAME=f8_chr13_brca2_gruhmb_age.vcf
SAVE_FILES=True
GNOMAD_FILE=gnomad_genome_31_sites_chr13_brca2.vcf

TOP_DIR=/private/groups/patenlab/jcasalet/federated-analysis
PYTHON_SCRIPT=${TOP_DIR}/app/cooccurrence/cooccurrenceFinder.py

miniwdl cromwell /private/groups/patenlab/jcasalet/WDL/myVusCooccur.wdl \
PYTHON_SCRIPT=${PYTHON_SCRIPT} \
VCF_FILE=${DATA_DIR}/$VCF_FILENAME \
ANNO_FILE=${DATA_DIR}/freeze8_sample_annot_2020-07-07.txt \
VARIANT_PATHOGENICITY_FILE=${TOP_DIR}/data/brca-variants.tsv \
OUTPUT_FILENAME=${GENE}-cooccurrences.json \
ALL_FILENAME=${GENE}-all.json \
VPI_FILENAME=${GENE}-vpi.json \
IPV_FILENAME=${GENE}-ipv.json \
TOUT_FILENAME=${GENE}-tout.json \
HG_VERSION=38 \
ENSEMBL_RELEASE=99 \
PHASED=True \
P2=$P2 \
CHROM=$CHROM \
GENE=$GENE \
NUM_CORES=$(grep -c processor /proc/cpuinfo) \
SAVE_FILES=$SAVE_FILES \
GNOMAD_FILE=${TOP_DIR}/data/$GNOMAD_FILE \
-c /private/groups/patenlab/jcasalet/WDL/cromwell.local.conf 
