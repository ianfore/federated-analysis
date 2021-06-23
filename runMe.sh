#!/bin/bash

FDA_PATH=$(pwd)
APP_PATH=${FDA_PATH}/app
CONF_PATH=/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/config
DATA_PATH=/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis/data
COOCCUR_DOCKER_IMAGE_NAME=brcachallenge/federated-analysis:cooccurrence
PATHOLOGY_DOCKER_IMAGE_NAME=brcachallenge/federated-analysis:pathology

if [ $# -eq 0 ]
then
	echo "usage: $0 -rc reportConfigFile -vf vcfFile -hg hgVersion -er ensemblRelease -c chromosome -p phased -g gene -pf pathogenicityFile -dd dataDirectory -st saveTempfiles -sp samplePathologyFile"
	echo "example: $0 -dq quality-report-config.json -vf my.vcf -hg 38 -er 99 -c 13 -p True --p2 0.001 -g BRCA2 -pf brca-variants.tsv -dd data -st True -sp mypf.tsv " 
	exit 1
	

else
	for arg in "$@"
	do
    		case $arg in
			-rc |--reportConfig)
			CONFIG_FILE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
			docker build -t ${PATHOLOGY_DOCKER_IMAGE_NAME} - < docker/pathology/Dockerfile
			docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}/pathology:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${PATHOLOGY_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py -c /config/$CONFIG_FILE 
			;;
			

			-vf|--vcfFile)
			VCF_FILE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-hg|--hgVersion)
			HG_VERSION="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-er|--ensemblRelease)
			ENSEMBL_RELEASE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-c|--chromosome)
			CHROM="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

        		-p|--phased)
        		PHASED="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

        		-p2|--p2)
        		P2="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-g|--gene)
			GENE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-pf|--pathogenicityFile)
			PATHOGENICITY_FILE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;


			-dd|--dataDirectory)
			DATA_PATH="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-st|--saveTempfiles)
			SAVE_FILES="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-sp|--samplePathologyFile)
			PATHOLOGY_FILE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

        		*)
			#echo "wrong usage: $0 $@"
			#exit 1
        		;;
    		esac
	done	


	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder.py  --vcf $VCF_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c $CHROM --p $PHASED --p2 $P2  --g $GENE --b $PATHOGENICITY_FILE  --d /var/tmp/pyensembl-cache  --data /data --save $SAVE_FILES --pf $PATHOLOGY_FILE


fi
