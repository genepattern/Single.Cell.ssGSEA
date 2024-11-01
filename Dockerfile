### copyright 203-2023. GenePattern Team @ Mesirov Lab - University of California, San Diego. All rights reserved.

FROM ubuntu:24.04

SHELL ["/bin/bash", "-c"]

WORKDIR /home/ubuntu/

## Prep for package installation
RUN apt update
RUN apt -y upgrade
RUN apt -y install sudo

## Configure users
RUN echo "ubuntu ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/ubuntu

## Install Python 3.12
RUN apt -y install build-essential zlib1g-dev \
libncurses5-dev libgdbm-dev libnss3-dev libssl-dev \
libreadline-dev libffi-dev libsqlite3-dev wget libbz2-dev \
pkg-config libhdf5-dev libpng-dev libxml2-dev liblapack-dev \
libopenblas-dev gfortran curl

RUN wget https://www.python.org/ftp/python/3.11.9/Python-3.11.9.tgz
RUN tar -xf Python-3.11.9.tgz
WORKDIR /home/ubuntu/Python-3.11.9
RUN ./configure --enable-optimizations
RUN make -j$(nproc)
RUN make install
WORKDIR /home/ubuntu/
RUN rm -r Python-3.11.9.tgz Python-3.11.9/

## Install pip
RUN apt -y install python3-pip

## Install Python3.11 venv
#RUN apt -y install python3.12-venv
RUN python3 -m pip install --upgrade virtualenv

## Install sc_ssGSEA
USER ubuntu
RUN python3 -m virtualenv sc_ssgsea_venv
ARG CACHEBUST=1 
ARG REPOSITORY_URL=https://pypi.org/simple
ARG VERSION=1.0.7
RUN source sc_ssgsea_venv/bin/activate && \
	python3 -m pip install --upgrade setuptools
RUN source sc_ssgsea_venv/bin/activate && \
	python3 -m pip install --index-url $REPOSITORY_URL --extra-index-url https://pypi.org/simple sc_ssGSEA==$VERSION

## Add driver from repo
COPY run_sc_ssgsea.py /home/ubuntu/
USER root
RUN chmod 777 /home/ubuntu/run_sc_ssgsea.py
USER ubuntu
