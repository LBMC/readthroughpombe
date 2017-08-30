#!/usr/bin/python3
# -*-coding:Utf-8 -*
# Copyright laurent modolo for the LBMC UMR 5239 ©.
# contributor(s) : laurent modolo (2017)
#
# laurent.modolo@ens-lyon.fr
#
# This software is a computer program whose purpose is to manage bioinformatic
# process in the pipeline of this project
#
# This software is governed by the CeCILL  license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

import sys
import argparse
import numpy
if sys.version_info[0] == 2:
    print("file_handle.py is only compatible with python3.\
    Please run file_handle.py as an executable or with the command\
    'python3 file_handle.py'")
    exit(1)

parser = argparse.ArgumentParser(
    prog="file_handle.py",
    description="script to transfrom kallisto pseudo outputs files into a \
    genes counts tables.")
parser.add_argument(
    "-t", "--transcripts",
    help="transcripts file used for the index (in fasta format)",
    type=str,
    action="store",
    dest="transcripts_file_name",
    required=True)
parser.add_argument(
    "-o", "--ouput_dir",
    help="kallisto pseudo output directory",
    type=str,
    action="store",
    dest="kallisto_output",
    required=True)
parser.add_argument(
    "-v",
    help="version information",
    default=False,
    action="store_true",
    dest="version")

args = parser.parse_args()

if args.version:
    print("0.0.1")
    exit(0)

cells_id = list()
with open(str(args.kallisto_output)+'/matrix.cells', 'r') as cells_file:
    for line in cells_file:
        cells_id.append(line[:-1])

transcripts = list()
with open(str(args.transcripts_file_name), 'r') as transcripts_file:
    for line in transcripts_file:
        if line[0] == '>':
            line = line.split("|")
            gene_name = line[5]
            transcripts.append(gene_name)

genes = list()
genes_pos = dict()
for transcript in transcripts:
    if transcript not in genes_pos:
        genes.append(transcript)
        genes_pos[transcript] = len(genes) - 1

transcripts_class = dict()
with open(str(args.kallisto_output)+'/matrix.ec', 'r') as transcripts_file:
    for line in transcripts_file:
        line = line.split()
        tclass = int(line[0])
        tlist = line[1].split(",")
        tlist = [int(n) for n in tlist]
        transcripts_class[tclass] = tlist

counts = numpy.zeros((len(cells_id), len(transcripts)))
with open(str(args.kallisto_output)+'/matrix.tsv', 'r') as counts_file:
    for line in counts_file:
        line = line.split()
        line = [int(n) for n in line]
        cell_pos = int(line[1])
        transcript_pos = int(line[0])
        transcript_class_counts = float(line[2])
        for transcript in transcripts_class[transcript_pos]:
            counts[cell_pos, transcript] += \
                transcript_class_counts / \
                float(len(transcripts_class[transcript_pos]))

with open(str(args.kallisto_output)+'/transcripts.counts', 'w') as output_file:
    output_file.write("id")
    for cell in range(0, len(cells_id)-1):
        output_file.write("\t"+cells_id[cell])
    output_file.write("\n")
    for transcript in range(0, len(transcripts)-1):
        output_file.write(transcripts[transcript])
        for cell in range(0, len(cells_id)-1):
            output_file.write("\t"+str(counts[cell, transcript]))
        output_file.write("\n")

counts = numpy.zeros((len(cells_id), len(genes)))
with open(str(args.kallisto_output)+'/matrix.tsv', 'r') as counts_file:
    for line in counts_file:
        line = line.split()
        line = [int(n) for n in line]
        cell_pos = int(line[1])
        transcript_pos = int(line[0])
        transcript_class_counts = float(line[2])
        for transcript in transcripts_class[transcript_pos]:
            gene_pos = genes_pos[transcripts[transcript]]
            counts[cell_pos, gene_pos] += \
                transcript_class_counts / \
                float(len(transcripts_class[transcript_pos]))

with open(str(args.kallisto_output)+'/genes.counts', 'w') as output_file:
    output_file.write("id")
    for cell in range(0, len(cells_id)-1):
        output_file.write("\t"+cells_id[cell])
    output_file.write("\n")
    for gene in range(0, len(genes)-1):
        output_file.write(genes[gene])
        for cell in range(0, len(cells_id)-1):
            output_file.write("\t"+str(counts[cell, gene]))
        output_file.write("\n")
