#!/bin/bash

FDA_PATH=/Users/jcasaletto/PycharmProjects/BIOBANK/federated-analysis

APP_PATH=${FDA_PATH}/app
CONF_PATH=${FDA_PATH}/config
DATA_PATH=${FDA_PATH}/data

#DOCKER_IMAGE_NAME=my_fda_image
DOCKER_IMAGE_NAME=mytest

docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py /config/conf.json

