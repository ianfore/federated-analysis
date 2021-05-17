#!/bin/bash


DATA_DIR=/private/groups/patenlab/jcasalet/TOPMED_WORKDIR/DATA/PUB/MANDY/BROAD/13/AMR
CHROM=13
GENE=BRCA2
P2=0.001
VCF_FILENAME=f8_chr13_brca2_mc_broad_amr.vcf
SAVE_FILES=True

TOP_DIR=/private/groups/patenlab/jcasalet/federated-analysis
PYTHON_SCRIPT=${TOP_DIR}/app/cooccurrence/cooccurrenceFinder.py
miniwdl cromwell /private/groups/patenlab/jcasalet/WDL/myVusCooccur.wdl \
PYTHON_SCRIPT=${PYTHON_SCRIPT} \
VCF_FILE=${DATA_DIR}/$VCF_FILENAME \
ANNO_FILE=${DATA_DIR}/freeze8_sample_annot_2020-07-07.txt \
BRCA_FILE=${TOP_DIR}/data/brca-variants.tsv \
OUTPUT_FILENAME=${CHROM}-out.json \
ALL_FILENAME=${CHROM}-all.json \
VPI_FILENAME=${CHROM}-vpi.json \
IPV_FILENAME=${CHROM}-ipv.json \
TOUT_FILENAME=${CHROM}-tout.json \
HG_VERSION=38 \
ENSEMBL_RELEASE=99 \
PHASED=True \
P2=$P2 \
CHROM=$CHROM \
GENE=$GENE \
NUM_CORES=$(grep -c processor /proc/cpuinfo) \
SAVE_FILES=$SAVE_FILES \
-c /private/groups/patenlab/jcasalet/WDL/cromwell.local.conf 
