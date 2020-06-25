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
	echo "$0 <input-vcf-filename> <output-json-filename> <hg-version> <ensembl-version> <chromosome-of-interest> <phased-boolean> <gene-of-interest> <save-variants> <brca-vars-filename>" 
	echo "example: $0 /data/bc3.vcf /data/myout.json 38 99 13 True BRCA2 True /data/brca-variants.tsv" 
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
	SAVE_VARS=$8
        BRCA_VARS=$9
	IPV_FILE=${10}
	VPI_FILE=${11}


	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder.py  --v ${VCF_FILE} --o ${OUTPUT_FILE} --h $HG_VERSION --e $ENSEMBL_RELEASE --c $CHROM --p $PHASED --s ${SAVE_VARS} --g $GENE --b $BRCA_VARS  --d /var/tmp/pyensembl-cache --ipv $IPV_FILE  --vpi $VPI_FILE



elif [ $# -eq 2 -a $1 == 'interactive' ]
then
	if [ $2 == 'root' ]
	then
		docker run -it --rm --user=0:0 -w /home/myuser -v "$(pwd)":/app  -v "${DATA_PATH}":/data:rw --entrypoint /bin/bash  ${COOCCUR_DOCKER_IMAGE_NAME}

	else
		docker run -it --rm --user=1968:1968 -w /home/myuser -v "$(pwd)":/app  -v "${DATA_PATH}":/data:rw --entrypoint /bin/bash  ${COOCCUR_DOCKER_IMAGE_NAME}
	fi

else
	echo "wrong usage"
	exit 1
	
fi
