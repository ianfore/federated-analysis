FROM ubuntu:18.04

MAINTAINER James Casaletto <james.casaletto@ucsc.edu>

USER root

RUN chmod 1777 /tmp /var/tmp

RUN useradd -u 1968 -g games -s /bin/bash -d /home/myuser -m myuser 

RUN apt-get update

RUN apt-get install -y \
    python3.7 \
    python3-pip \
    python3-dev 


#RUN pip3 install --no-cache-dir numpy pandas pandasql sklearn tabulate scipy pyensembl scikit-allel hmmlearn matplotlib
RUN pip3 install --no-cache-dir numpy pandas pyensembl scikit-allel 


RUN pyensembl install --release 75 --species homo_sapiens
RUN pyensembl install --release 99 --species homo_sapiens

RUN mv /root/.cache /var/tmp/pyensembl-cache
RUN chmod 755 /var/tmp/pyensembl-cache
