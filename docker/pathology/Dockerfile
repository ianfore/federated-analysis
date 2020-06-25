FROM ubuntu:18.04

MAINTAINER James Casaletto <james.casaletto@ucsc.edu>

USER root

RUN chmod 1777 /tmp /var/tmp

RUN useradd -u 1968 -g games -s /bin/bash -d /home/myuser -m myuser 

RUN apt-get update && apt-get install -y  python3.7  python3-pip  python3-dev software-properties-common && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && DEBIAN_FRONTEND=noninteractive apt-get install r-base -y

RUN pip3 install --no-cache-dir numpy pandas pandasql sklearn tabulate scipy pyensembl scikit-allel hmmlearn




RUN pyensembl install --release 75 --species homo_sapiens
RUN pyensembl install --release 99 --species homo_sapiens

RUN mv /root/.cache /var/tmp/pyensembl-cache
RUN chmod 755 /var/tmp/pyensembl-cache
