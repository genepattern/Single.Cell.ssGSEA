### copyright 203-2023. GenePattern Team @ Mesirov Lab - University of California, San Diego. All rights reserved.

FROM satijalab/seurat:5.0.0

SHELL ["/bin/bash", "-c"]

## Install Python 3.11.9

RUN wget https://www.python.org/ftp/python/3.11.9/Python-3.11.9.tgz
RUN tar -xf Python-3.11.9.tgz
WORKDIR /Python-3.11.9
RUN ./configure --enable-optimizations
RUN make -j$(nproc)
RUN make install
WORKDIR /
RUN rm -r Python-3.11.9.tgz Python-3.11.9/

## Install pip

RUN apt -y install python3-pip

## Install Python3.11 venv

RUN python3 -m pip install --upgrade virtualenv

## Install sc_ssGSEA

ARG CACHEBUST=1 
RUN python3 -m virtualenv sc_ssgsea_venv
RUN source sc_ssgsea_venv/bin/activate && \
	python3 -m pip install --upgrade setuptools
COPY ./ sc_ssGSEA_local
RUN source sc_ssgsea_venv/bin/activate && \
	python3 -m pip install /sc_ssGSEA_local

## Add driver from repo
COPY run_sc_ssgsea.py /
RUN chmod 777 run_sc_ssgsea.py

CMD ["/bin/bash && source sc_ssgsea_venv/bin/activate"]
