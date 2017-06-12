FROM ubuntu:16.04
MAINTAINER Laurent Modolo

RUN apt-get update
RUN apt-get install -y wget git cmake python python3-pip trimmomatic
RUN pip3 install argparse cutadapt
RUN pip install multiqc click jinja2 htseq

#install UrQt
RUN git clone git@github.com:l-modolo/UrQt.git UrQt
RUN cd UrQt
RUN make
RUN cp UrQt /usr/local/bin/
