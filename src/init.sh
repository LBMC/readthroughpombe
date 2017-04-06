#!/bin/sh

# install nextflow
cd bin
wget -qO- get.nextflow.io | bash
pip install docker
pip install multiqc
