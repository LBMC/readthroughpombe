FROM ubuntu:16.04
MAINTAINER Laurent Modolo

RUN apt-get update && apt-get install -y \
    wget=1.17.1-1ubuntu1.2 \
    git=1:2.7.4-0ubuntu1.1 \
    python-pip=8.1.1-2ubuntu0.4 \
    python3-pip=8.1.1-2ubuntu0.4 \
    build-essential=12.1ubuntu2 \

   && rm -rf /var/lib/apt/lists/*
RUN pip3 install argparse==1.4.0
RUN pip install numpy==1.12.1 && \
  pip install multiqc==0.9 htseq==0.7.2

RUN git clone https://github.com/COMBINE-lab/salmon.git && \
  cd salmon && \
  git checkout tags/v0.8.2 -b master && \
  mkdir build && \
  cmake -DFETCH_BOOST=TRUE -DCMAKE_INSTALL_PREFIX=/usr/  && \
  make  && \
  sudo make install && \
  rm -Rf /tmp/salmon
