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
	echo "$0 <input-vcf-filename> <hg-version> <ensembl-version> <chromosome-of-interest> <phased-boolean> <gene-of-interest> <brca-vars-filename> <save-files-boolean>" 
	echo "example: $0 bc3.vcf 38 99 13 True BRCA2 brca-variants.tsv True" 
	exit 1
	

elif [ $# -eq 1 -a $1 == 'analyze' ]
then

	docker build -t ${PATHOLOGY_DOCKER_IMAGE_NAME} - < docker/pathology/Dockerfile

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}/pathology:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${PATHOLOGY_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py /config/conf.json 
	exit 0

elif [ $# -eq 2 -a $1 == 'interactive' ]
then
	if [ $2 == 'root' ]
	then
		docker run -it --rm --user=0:0 -w /home/myuser -v "$(pwd)":/app  -v "${DATA_PATH}":/data:rw --entrypoint /bin/bash  ${COOCCUR_DOCKER_IMAGE_NAME}
	else
		docker run -it --rm --user=1968:1968 -w /home/myuser -v "$(pwd)":/app  -v "${DATA_PATH}":/data:rw --entrypoint /bin/bash  ${COOCCUR_DOCKER_IMAGE_NAME}
	fi

elif [ $# -eq 9 ]
then
	
	VCF_FILE=$1
	HG_VERSION=$2
	ENSEMBL_RELEASE=$3
	CHROM=$4
	PHASED=$5
	GENE=$6
        BRCA_VARS=$7
	DATA_DIR=$8
	SAVE_FILES=$9


	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder.py  --vcf $VCF_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c $CHROM --p $PHASED  --g $GENE --b $BRCA_VARS  --d /var/tmp/pyensembl-cache  --data /data --save $SAVE_FILES


else
	echo "wrong usage"
	exit 1
fi
