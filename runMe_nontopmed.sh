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
	echo "$0 <input-vcf-filename> <chromosome-of-interest> <gene-of-interest>  controlsOnly" 
	echo "example: $0 breastcancer.vcf 13 BRCA2 controlsOnly" 
	echo "OR"
	echo "$0 <input-vcf-filename> <chromosome-of-interest> <gene-of-interest>  casesOnly <pathology-report>" 
	echo "example: $0 breastcancer.vcf 13 BRCA2 casesOnly pathology.tsv" 
	exit 1
	

elif [ $# -eq 1 -a $1 == 'analyze' ]
then
	docker build -t ${PATHOLOGY_DOCKER_IMAGE_NAME} - < docker/pathology/Dockerfile
	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}/pathology:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${PATHOLOGY_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py /config/conf.json 
	exit 0

elif [ $# -eq 4 ]
then
	VCF_FILE=/data/$1
	CHROM=$2
	GENE=$3
	CASES_OR_CONTROLS=$4
	HG_VERSION=37
	ENSEMBL_RELEASE=75
	PHASED=False
        BRCA_VARS=/data/brca-variants.tsv
	ALL_FILE=/data/${CHROM}-all-controlsOnly.json
	OUTPUT_FILE=/data/${CHROM}-out-controlsOnly.json

	if [ "$CASES_OR_CONTROLS" == "casesOnly" ]
	then
		echo "you must provide a pathology file when running casesOnly"
		exit 1
	elif [ "$CASES_OR_CONTROLS" != "controlsOnly" ]
	then
		echo "must use controlsOnly when not providing a pathology file"
		exit 1
	fi

	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder_nontopmed.py --vcf $VCF_FILE --out $OUTPUT_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c $CHROM --g $GENE --p $PHASED --b $BRCA_VARS --all $ALL_FILE --controlsOnly


elif [ $# -eq 5 ]
then
	VCF_FILE=/data/$1
	CHROM=$2
	GENE=$3
	CASES_OR_CONTROLS=$4
	PATHOLOGY_FILE=/data/$5
	HG_VERSION=37
	ENSEMBL_RELEASE=75
	PHASED=False
        BRCA_VARS=/data/brca-variants.tsv
	ALL_FILE=/data/${CHROM}-all-casesOnly.json
	OUTPUT_FILE=/data/${CHROM}-out-casesOnly.json

	if [ "$CASES_OR_CONTROLS" == "controlsOnly" ]
	then
		echo "do not provide a pathology file when running controlsOnly"
		exit 1
	elif [ "$CASES_OR_CONTROLS" != "casesOnly" ]
	then
		echo "must use casesOnly when providing a pathology file"
		exit 1
	fi

	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder_nontopmed.py --vcf $VCF_FILE --out $OUTPUT_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c $CHROM --g $GENE --p $PHASED --b $BRCA_VARS --all $ALL_FILE --pf $PATHOLOGY_FILE --casesOnly

else
	echo "wrong usage"
	exit 1
fi
