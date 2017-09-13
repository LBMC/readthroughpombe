#!/bin/sh

# This script renames files using file_handle. 
# bash date.sh "file"

python3 ./src/file_handle/src/file_handle.py -f $1
