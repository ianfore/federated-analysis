#!/bin/bash

DATA_DIR=/private/groups/patenlab/jcasalet/TOPMED_WORKDIR/DATA/BROAD/13/GNOMAD/MC
TOP_DIR=/private/groups/patenlab/jcasalet/federated-analysis
PYTHON_SCRIPT=${TOP_DIR}/app/cooccurrence/cooccurrenceFinder.py

miniwdl cromwell /private/groups/patenlab/jcasalet/WDL/myVusCooccur.wdl \
PYTHON_SCRIPT=${PYTHON_SCRIPT} \
VCF_FILE=${DATA_DIR}/NOTINSUBSET/broad_chr13_brca2_notinsubset.vcf \
BRCA_FILE=${TOP_DIR}/data/brca-variants.tsv \
OUTPUT_FILENAME=13-notinsubset-out.json \
HG_VERSION=38 \
ENSEMBL_RELEASE=99 \
PHASED=True \
CHROM=13 \
GENE=BRCA2 \
NUM_CORES=$(grep -c processor /proc/cpuinfo) \
VPI_FILENAME=13-notinsubset-vpi.json \
IPV_FILENAME=13-notinsubset-ipv.json \
-c /private/groups/patenlab/jcasalet/WDL/cromwell.local.conf 
