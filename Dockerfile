FROM ubuntu:18.04

MAINTAINER James Casaletto <james.casaletto@ucsc.edu>

USER root


RUN apt-get update && apt-get install -y \
    python3.5 \
    python3-pip \
    python3-dev 


RUN pip3 install --no-cache-dir numpy pandas pandasql tabulate scipy

RUN apt-get install software-properties-common -y

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

RUN DEBIAN_FRONTEND=noninteractive  apt-get install r-base -y
