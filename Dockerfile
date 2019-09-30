FROM ubuntu:18.04

MAINTAINER James Casaletto <james.casaletto@ucsc.edu>

USER root

RUN apt-get update && apt-get install -y \
    python3.5 \
    python3-pip \
    python3-dev 


RUN pip3 install --no-cache-dir numpy pandas pandasql tabulate 

