#!/bin/bash

DATA_DIR=/private/groups/patenlab/jcasalet/TOPMED_WORKDIR/DATA/BROAD/13/GNOMAD/MC
TOP_DIR=/private/groups/patenlab/jcasalet/federated-analysis
PYTHON_SCRIPT=${TOP_DIR}/app/cooccurrence/cooccurrenceFinder.py

miniwdl cromwell /private/groups/patenlab/jcasalet/WDL/myVusCooccur.wdl \
PYTHON_SCRIPT=${PYTHON_SCRIPT} \
VCF_FILE=${DATA_DIR}/NOTINSUBSET/broad_chr13_brca2_notinsubset.vcf \
BRCA_FILE=${TOP_DIR}/data/brca-variants.tsv \
OUTPUT_FILENAME=13-out.json \
ALL_FILENAME=13-all.json \
VPI_FILENAME=13-vpi.json \
IPV_FILENAME=13-ipv.json \
HG_VERSION=38 \
FREEZE=8 \
ENSEMBL_RELEASE=99 \
PHASED=True \
CHROM=13 \
GENE=BRCA2 \
NUM_CORES=$(grep -c processor /proc/cpuinfo) \
-c /private/groups/patenlab/jcasalet/WDL/cromwell.local.conf 
