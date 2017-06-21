FROM ubuntu:16.04
MAINTAINER Laurent Modolo

ENV SALMON_VERSION 0.8.2
ENV KALLISTO_VERSION 0.43.1
ENV BOWTIE2_version 2.3.2

ENV PACKAGES curl=7.47.0-1ubuntu2.2 \
    git=1:2.7.4-0ubuntu1.1 \
    unzip=6.0-20ubuntu1 \
    pigz=2.3.1-2 \
    ca-certificates=20160104ubuntu1 \
    python-pip=8.1.1-2ubuntu0.4 \
    python3-pip=8.1.1-2ubuntu0.4 \
    g++=4:5.3.1-1ubuntu1 \
    gcc=4:5.3.1-1ubuntu1 \
    autoconf=2.69-9 \
    pkg-config=0.29.1-0ubuntu1 \
    make=4.1-6 \
    cmake=3.5.1-1ubuntu3 \
    libboost-all-dev=1.58.0.1ubuntu1 \
    liblzma-dev=5.1.1alpha+20120614-2ubuntu2 \
    libbz2-dev=1.0.6-8 \
    zlib1g-dev=1:1.2.8.dfsg-2ubuntu4.1 \
    libjemalloc-dev=3.6.0-9ubuntu1 \
    libdivsufsort-dev=2.0.1-2

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN pip3 install argparse==1.4.0
RUN pip install setuptools==36.0.1 && \
    pip install multiqc==0.9


RUN  git clone https://github.com/pachterlab/kallisto.git && \
  cd kallisto && \
  git checkout tags/v0.43.1 && \
  mkdir build && \
  cmake && \
  make && \
  make install && \
  cd .. && rm -Rf kallisto

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip/download && \
  unzip download && \
  cd bowtie2-2.3.2 && \
  cp bowtie2* /usr/bin/ && \
  cd .. && rm -Rf bowtie2* download

RUN  git clone https://github.com/deweylab/RSEM.git && \
  cd RSEM && \
  git checkout tags/v1.3.0 && \
  make && \
  sudo make install && \
  cd .. && rm -Rf RSEM

COPY file_handle.py /usr/bin/local/file_handle.py
