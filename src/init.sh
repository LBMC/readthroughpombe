#!/bin/sh

# install nextflow
cd bin
wget -qO- get.nextflow.io | bash
pip install argparse
pip install docker
pip2 install multiqc
pip2 install click
pip2 install jinja2
