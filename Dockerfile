##### BASE IMAGE #####
FROM python:3.8-slim

##### METADATA #####
LABEL base.image="python:3.8-slim"
LABEL software="MultiQC + custom plugins"
LABEL software.description="Aggregate bioinformatics results in single report"
LABEL software.website="https://github.com/zavolanlab/multiqc-plugins"
LABEL software.documentation="https://github.com/zavolanlab/multiqc-plugins"
LABEL software.license="https://github.com/zavolanlab/multiqc-plugins/blob/master/LICENSE"
LABEL software.tags="Bioinformatcs"
LABEL maintainer="maciej.bak@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="Apache License 2.0"

COPY modules modules
COPY tests tests
COPY setup.py .

##### INSTALL #####
RUN apt-get update \
  && apt-get install -y apt-utils gcc make zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
  && pip install multiqc==1.9 \
  && python setup.py install
