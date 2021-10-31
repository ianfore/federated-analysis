#!/bin/bash

FDA_PATH=$(pwd)
APP_PATH=${FDA_PATH}/app
COOCCUR_DOCKER_IMAGE_NAME=brcachallenge/federated-analysis:cooccurrence
PATHOLOGY_DOCKER_IMAGE_NAME=brcachallenge/federated-analysis:pathology

if [ $# -eq 0 ]
then
	echo "example: $0 -c 13 -p True -dd $(pwd)/examples/BRCA2/data -cd $(pwd)/examples/BRCA2/config -vf brca2.vcf -vpf clinvar_brca2.tsv -gf gnomad_chr13_brca2.vcf -rc brca2-report-config.json -g BRCA2 -spf brca2-pathology.tsv -hg 38" 
	exit 1
	

else
	for arg in "$@"
	do
    		case $arg in
			-gf |--gnomadFile)
			GNOMAD_FILE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
			;;

			-rc |--reportConfig)
			CONFIG_FILE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
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

			-vpf|--variantPathogenicityFile)
			PATHOGENICITY_FILE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;


			-dd|--dataDirectory)
			DATA_PATH="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-cd|--confDirectory)
			CONF_PATH="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-st|--saveTempfiles)
			SAVE_FILES="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

			-spf|--samplePhenotypeFile)
			PHENOTYPE_FILE="$2"
        		shift # Remove argument name from processing
        		shift # Remove argument value from processing
        		;;

        		*)
			#echo "wrong usage: $0 $@"
			#exit 1
        		;;
    		esac
	done	

	if [ -n "$CONFIG_FILE" ]
	then
		docker build -t ${PATHOLOGY_DOCKER_IMAGE_NAME} - < docker/pathology/Dockerfile
		docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}/pathology:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${PATHOLOGY_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py -c /config/$CONFIG_FILE 
	fi

	docker build -t ${COOCCUR_DOCKER_IMAGE_NAME} - < docker/cooccurrence/Dockerfile

	if [ -z "$HG_VERSION" ]
	then
		HG_VERSION=38
	fi
	if [ -z "$ENSEMBL_RELEASE" ]
	then
		ENSEMBL_RELEASE=99
	fi
	if [ -z "$SAVE_FILES" ]
	then
		SAVE_FILES=False
	fi
	if [ -z "$P2" ]
	then
		P2=0.001	
	fi
	if [ -z "PHENOTYPE_FILE" ]
	then
		PHENOTYPE_FILE=""
	fi

	docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 --user=`id -u`:`id -g` -v ${APP_PATH}/cooccurrence:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${COOCCUR_DOCKER_IMAGE_NAME} /usr/bin/python3 /app/cooccurrenceFinder.py  --vcf $VCF_FILE --h $HG_VERSION --e $ENSEMBL_RELEASE --c $CHROM --p $PHASED --p2 $P2  --g $GENE --vpf $PATHOGENICITY_FILE  --d /var/tmp/pyensembl-cache  --data /data --save $SAVE_FILES --spf "$PHENOTYPE_FILE"  --gf "$GNOMAD_FILE"


fi
