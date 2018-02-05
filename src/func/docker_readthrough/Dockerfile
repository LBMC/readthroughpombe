FROM ubuntu:16.04
MAINTAINER Laurent Modolo

ENV FILE_HANDLE_VERSION 0.1.1
ENV BOWTIE2_VERSION 2.3.2
ENV SAMTOOLS_VERSION 1.5
ENV MUSIC_VERSION 6613c532bf21d958ae380f064cc7752f3e8c368f
ENV BEDOOPS_VERSION 2.4.30

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
    libfreetype6=2.6.1-0.1* \
    genometools=1.5* \
    bedtools=2.25.0-1

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

# install file_handle
RUN git clone https://github.com/LBMC/file_handle.git && \
    cd file_handle && \
    git checkout tags/v${FILE_HANDLE_VERSION} && \
    cp src/file_handle.py /usr/bin/file_handle.py && \
    chmod +x /usr/bin/file_handle.py

# install MUSIC
RUN git clone https://github.com/gersteinlab/MUSIC.git && \
    cd MUSIC && \
    git checkout ${MUSIC_VERSION} &&\
    make clean && \
    make && \
    chmod +x bin/* && \
    cp bin/* /usr/local/bin/

#install bowtie2
RUN curl -k -L https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip/download -o bowtie2-${BOWTIE2_VERSION}.zip && \
  unzip bowtie2-${BOWTIE2_VERSION}.zip && \
  cd bowtie2-${BOWTIE2_VERSION} && \
  cp bowtie2* /usr/local/bin/ && \
  cd .. && rm -Rf bowtie2*

# install samtools
RUN curl -k -L https://sourceforge.net/projects/samtools/files/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2/download -o samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
  tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
  cd samtools-${SAMTOOLS_VERSION} && \
  ./configure && \
  make && \
  make install && \
  cd .. && rm -Rf samtools-${SAMTOOLS_VERSION}*

# install bedops
RUN curl -k -L https://github.com/bedops/bedops/releases/download/v2.4.30/bedops_linux_x86_64-v${BEDOOPS_VERSION}.tar.bz2 -o bedops_linux_x86_64-v${BEDOOPS_VERSION}.tar.bz2 && \
  tar xvf bedops_linux_x86_64-v${BEDOOPS_VERSION}.tar.bz2 && \
  cp bin/* /usr/local/bin/

# install kallisto
ENV KALLISTO_VERSION 0.43.1
RUN  curl -k -L https://github.com/pachterlab/kallisto/releases/download/v${KALLISTO_VERSION}/kallisto_linux-v${KALLISTO_VERSION}.tar.gz -o kallisto_linux-v${KALLISTO_VERSION}.tar.gz && \
  tar xzf kallisto_linux-v${KALLISTO_VERSION}.tar.gz && \
  cp kallisto_linux-v${KALLISTO_VERSION}/kallisto /usr/local/bin/ && \
  cd .. && rm -Rf kallisto_linux-v${KALLISTO_VERSION}*
