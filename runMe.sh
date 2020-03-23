#!/bin/bash

FDA_PATH=$(pwd)
APP_PATH=${FDA_PATH}/app
CONF_PATH=${FDA_PATH}/config
DATA_PATH=${FDA_PATH}/data
DOCKER_IMAGE_NAME=brcachallenge/federated-analysis

if [ $# -eq 0 ]
then
	echo "usage: $0 <vcf-file> analyze"
	echo "OR"
	echo "$0 <input-vcf-filename> <output-json-filename> <hg-version> <ensembl-version> <chromosomes-of-interest> <phased-boolean>"
	echo "example: $0 BreastCancer.shuffle.vcf myCooccurrences.json 37 75 False"
	exit 1
	

elif [ $# -eq 1 -a $1 == 'analyze' ]
then
	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py /config/conf.json 
	exit 0
elif [ $# -eq 5 ]
then

	VCF_FILE=$1
	OUTPUT_FILE=$2
	HG_VERSION=$3
	ENSEMBL_RELEASE=$4
	PHASED=$5

	CHR_LIST=[13,17]

	if [ ! -d ${DATA_PATH}/pyensembl-cache ]
	then
        	cd $DATA_PATH
        	cat pyens.?? > pyensembl.zip
        	unzip pyensembl.zip
		rm pyens.??
		rm pyensembl.zip
        	cd -
	fi

	if [ $(uname) == "Darwin" ]
	then
		PREV_PERMS=$(stat -f "%OLp" ${DATA_PATH})
	else
		PREV_PERMS=$(stat -c "%a" ${DATA_PATH})
	fi

	chmod 1777 ${DATA_PATH}

	docker build -t ${DOCKER_IMAGE_NAME} .

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w /home/myuser --user=1968:games -v ${APP_PATH}:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder.py  ${VCF_FILE} ${OUTPUT_FILE} --h $HG_VERSION --e $ENSEMBL_RELEASE --c $CHR_LIST --p $PHASED

	#docker run -it --rm --user=1968:1968 -w /home/myuser -v "$(pwd)":/app  -v "${DATA_PATH}":/data:rw --entrypoint /bin/bash  ${DOCKER_IMAGE_NAME}
	#docker run -it --rm --user=0:0 -w /home/myuser -v "$(pwd)":/app  -v "${DATA_PATH}":/data:rw --entrypoint /bin/bash  ${DOCKER_IMAGE_NAME}

	chmod $PREV_PERMS ${DATA_PATH}
fi
