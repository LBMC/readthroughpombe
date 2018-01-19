FROM ubuntu:16.04
MAINTAINER Laurent Modolo

ENV FILE_HANDLE_VERSION 0.1.1
ENV SEQKIT_VERSION 0.7.2

ENV PACKAGES curl=7.47.0* \
    git=1:2.7.4-0* \
    unzip=6.0-20* \
    pigz=2.3.1-2 \
    ca-certificates=20160104* \
    python-pip=8.1.1-2* \
    python3-pip=8.1.1-2* \
    g++=4:5.3.1-1* \
    gcc=4:5.3.1-1* \
    autoconf=2.69-9 \
    pkg-config=0.29.1-0* \
    make=4.1-6 \
    cmake=3.5.1-1* \
    libboost-all-dev=1.58.0.1* \
    liblzma-dev=5.1.1alpha+20120614-2* \
    libbz2-dev=1.0.6-8 \
    zlib1g-dev=1:1.2.8.dfsg-2* \
    libjemalloc-dev=3.6.0-9* \
    libdivsufsort-dev=2.0.1-2 \
    libtbb2=4.4~20151115-0* \
    libncurses5=6.0+20160213-1* \
    libncurses5-dev=6.0+20160213-1* \
    libfreetype6=2.6.1-0.1*

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -o seqkit_linux_amd64.tar.gz && \
  tar xvf seqkit_linux_amd64.tar.gz && \
  cp seqkit /usr/local/bin/

# install file_handle
RUN git clone https://github.com/LBMC/file_handle.git && \
    cd file_handle && \
    git checkout tags/v${FILE_HANDLE_VERSION} && \
    cp src/file_handle.py /usr/bin/file_handle.py && \
    chmod +x /usr/bin/file_handle.py
