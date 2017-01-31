# Build as jupyter/systemuser
# Run with the DockerSpawner in JupyterHub

FROM continuumio/anaconda

MAINTAINER ENCODE-ELTE-NODE

USER root
ENV SHELL /bin/bash

#########################################
#baseline dependecies
#########################################
#pybedtools pybigwig
RUN conda install -y -c bioconda pybedtools pybigwig
RUN conda install -y psycopg2=2.6.2 biopython shapely
RUN apt-get update && apt-get install -y bedtools gcc

#pydnabinding
RUN git clone https://github.com/nboley/pyDNAbinding.git
RUN cd pyDNAbinding; python setup.py build; python setup.py install
#install failed to copy this directory
RUN cp -r /pyDNAbinding/pyDNAbinding/shape_data/ /opt/conda/lib/python2.7/site-packages/pyDNAbinding-0.1.1-py2.7-linux-x86_64.egg/pyDNAbinding/

#########################################
# other dependecies
#########################################
RUN apt-get install -y sudo screen
RUN apt-get install samtools

#pydnase
RUN pip install pyDNase
#Millipede ? 

#########################################
# jupyterhub dependecies
#########################################
#create 3.4 conda env
RUN conda create -n py34 python=3.4
RUN /bin/bash -c "source activate py34 && \
	conda install -y -c conda-forge jupyterhub && \
	conda install -y notebook "
#add python2 kernel to jupyter
RUN /bin/bash -c "source activate py34 && \
	python2 -m ipykernel install "

#########################################
# R and dependencies
#########################################
# install R
RUN conda install -c r r-essentials
# update irkernel for newestone
RUN conda update -y -c r r-irkernel
# bug with R plotting in newest r-base, downgrade it to a working one
# https://github.com/jupyter/docker-stacks/issues/210#issuecomment-246081809
RUN conda install -y -c r r-base='3.3.1 1'


#########################################
# New stuff goes here for faster building
#########################################
