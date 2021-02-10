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
	echo "$0 -v <input-vcf-filename> -m controlsOnly" 
	echo "example: $0 -v breastcancer.vcf -m controlsOnly" 
	echo "OR"
	echo "$0 -v <input-vcf-filename> -m casesOnly -p <pathology-report>" 
	echo "example: $0 -v breastcancer.vcf -m casesOnly -p pathology.tsv" 
	exit 1
	

elif [ $# -eq 1 -a $1 == 'analyze' ]
then
	docker build -t ${PATHOLOGY_DOCKER_IMAGE_NAME} - < docker/pathology/Dockerfile
	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}/pathology:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${PATHOLOGY_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py /config/conf.json 
	exit 0

elif [ $# -eq 4 ]
then
	while getopts ":v:m:" opt
	do
		case $opt in
		v) VCF_FILE=/data/$OPTARG
		;;
		m) MODE=$OPTARG
		;;
		esac
	done
	echo "vcf = " $VCF_FILE
	echo "mode = " $MODE
	HG_VERSION=37
	ENSEMBL_RELEASE=75
	PHASED=False
        BRCA_VARS=/data/brca-variants.tsv

	if [ "$MODE" == "casesOnly" ]
	then
		echo "you must provide a pathology file when running casesOnly"
		exit 1
	elif [ "$MODE" != "controlsOnly" ]
	then
		echo "must use controlsOnly when not providing a pathology file"
		exit 1
	fi

	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	ALL_FILE=/data/13-all-controlsOnly.json
	OUTPUT_FILE=/data/13-out-controlsOnly.json
	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder_nontopmed.py --vcf $VCF_FILE --out $OUTPUT_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c 13 --g BRCA2 --p $PHASED --b $BRCA_VARS --all $ALL_FILE --controlsOnly

	ALL_FILE=/data/17-all-controlsOnly.json
	OUTPUT_FILE=/data/17-out-controlsOnly.json
	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder_nontopmed.py --vcf $VCF_FILE --out $OUTPUT_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c 17 --g BRCA1 --p $PHASED --b $BRCA_VARS --all $ALL_FILE --controlsOnly


elif [ $# -eq 6 ]
then
	while getopts ":v:m:p:" opt
	do
		case $opt in
		v) VCF_FILE=/data/$OPTARG
		;;
		m) MODE=$OPTARG
		;;
		p) PATHOLOGY_FILE=/data/$OPTARG
		;;
		esac
	done
	HG_VERSION=37
	ENSEMBL_RELEASE=75
	PHASED=False
        BRCA_VARS=/data/brca-variants.tsv

	if [ "$MODE" == "controlsOnly" ]
	then
		echo "do not provide a pathology file when running controlsOnly"
		exit 1
	elif [ "$MODE" != "casesOnly" ]
	then
		echo "must use casesOnly when providing a pathology file"
		exit 1
	fi

	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	ALL_FILE=/data/13-all-casesOnly.json
	OUTPUT_FILE=/data/13-out-casesOnly.json
	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder_nontopmed.py --vcf $VCF_FILE --out $OUTPUT_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c 13 --g BRCA2 --p $PHASED --b $BRCA_VARS --all $ALL_FILE --pf $PATHOLOGY_FILE --casesOnly

	ALL_FILE=/data/17-all-casesOnly.json
	OUTPUT_FILE=/data/17-out-casesOnly.json
	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder_nontopmed.py --vcf $VCF_FILE --out $OUTPUT_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c 17 --g BRCA1 --p $PHASED --b $BRCA_VARS --all $ALL_FILE --pf $PATHOLOGY_FILE --casesOnly

else
	echo "wrong usage"
	exit 1
fi
