#!/bin/bash

if [ $# -ne 1 -a $0 != "-p" -a $0 != "-c" ]
then
	echo "usage: $0 -p|-c|-b"
	exit 1
fi

OP=$1
FDA_PATH=$(pwd)
APP_PATH=${FDA_PATH}/app
CONF_PATH=${FDA_PATH}/config
DATA_PATH=${FDA_PATH}/data
DOCKER_IMAGE_NAME=brcachallenge/federated-analysis

PREV_PERMS=$(stat -c "%a" ${DATA_PATH})

chmod 1777 ${DATA_PATH}

docker build -t ${DOCKER_IMAGE_NAME} .

#docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=`id -u`:`id -g` -v ${APP_PATH}:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data ${DOCKER_IMAGE_NAME} /usr/bin/python3 /app/dataAnalyzer.py /config/conf.json 
docker run --rm -e PYTHONPATH=/ -e PYTHONIOENCODING=UTF-8 -w / --user=1968:games -v ${APP_PATH}:/app:ro -v ${CONF_PATH}:/config -v "${DATA_PATH}":/data:rw ${DOCKER_IMAGE_NAME} /usr/bin/python3 /app/myCooccurrenceFinder.py  $OP 

#docker run -it --rm --user=1968:1968 -v "$(pwd)":/app  -v "${DATA_PATH}":/data:rw --entrypoint /bin/bash  ${DOCKER_IMAGE_NAME}

chmod $PREV_PERMS ${DATA_PATH}
