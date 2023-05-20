### copyright 203-2023. GenePattern Team @ Mesirov Lab - University of California, San Diego. All rights reserved.
#
# Currently, module uses genepattern/seurat-suite:4.0.3 image.
FROM genepattern/seurat-suite:4.0.3

LABEL maintainer="John Jun johnjun094@cloud.ucsd.edu"

USER root

RUN pip3 install tqdm==4.65.0 numpy==1.23.1 matplotlib==3.6.1 scanpy==1.9.1 pandas==1.4.4 anndata==0.8.0 \
    seaborn==0.12.0 scipy==1.9.2 networkx==2.8.7 xlsxwriter==3.0.3 openpyxl==3.0.9 mygene==3.2.2 humanfriendly==10.0

RUN Rscript -e "install.packages('optparse', version='1.7.3', repos='http://cran.us.r-project.org')"

# COPY R and Python scripts
RUN mkdir /scripts
COPY run_scgsea.py preprocess.R scgsea_helper.py /scripts/
RUN chmod a+rwx /scripts/*
