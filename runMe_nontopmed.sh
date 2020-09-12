#!/bin/bash

FDA_PATH=$(pwd)
APP_PATH=${FDA_PATH}/app
CONF_PATH=${FDA_PATH}/config
DATA_PATH=${FDA_PATH}/data
COOCCUR_DOCKER_IMAGE_NAME=brcachallenge/federated-analysis:cooccurrence
PATHOLOGY_DOCKER_IMAGE_NAME=brcachallenge/federated-analysis:pathology

if [ $# -eq 0 ]
then
	echo "usage: $0 <vcf-file> analyze"
	echo "OR"
	echo "$0 <input-vcf-filename> <output-json-filename> <hg-version> <ensembl-version> <chromosome-of-interest> <phased-boolean> <gene-of-interest> <brca-vars-filename> <ipv-output-file> <vpi-output-file> <all-output-file>" 
	echo "example: $0 /data/breastcancer.vcf /data/13-out.json 37 75 13 False BRCA2 /data/brca-variants.tsv 13-ipv.json 13-vpi.json 13-all.json" 
	exit 1
	

elif [ $# -eq 1 -a $1 == 'analyze' ]
then

	docker build -t ${PATHOLOGY_DOCKER_IMAGE_NAME} - < docker/pathology/Dockerfile

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}/pathology:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${PATHOLOGY_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py /config/conf.json 
	exit 0

elif [ $# -eq 11 ]
then

	VCF_FILE=$1
	OUTPUT_FILE=$2
	HG_VERSION=$3
	ENSEMBL_RELEASE=$4
	CHROM=$5
	PHASED=$6
	GENE=$7
        BRCA_VARS=$8
	IPV_FILE=${9}
	VPI_FILE=${10}
	ALL_FILE=${11}


	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder_nontopmed.py --vcf $VCF_FILE --out $OUTPUT_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c $CHROM --g $GENE --p $PHASED --n 1 --b $BRCA_VARS --vpi $VPI_FILE --ipv $IPV_FILE --all $ALL_FILE 


else
	echo "wrong usage"
	exit 1
fi
