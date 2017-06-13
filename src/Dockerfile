FROM ubuntu:14.04
MAINTAINER Laurent Modolo

RUN apt-get update
RUN apt-get install -y wget git cmake python-pip build-essential python3-pip trimmomatic zlib1g-dev libbz2-dev liblzma-dev
RUN pip3 install --upgrade pip
RUN pip3 install argparse encoding
RUN pip install numpy click jinja2
RUN pip install multiqc htseq
RUN export LC_ALL=en_US.UTF8 && pip install cutadapt

#install UrQt
RUN git clone git@github.com:l-modolo/UrQt.git UrQt
RUN cd UrQt
RUN make
RUN cp UrQt /usr/local/bin/
