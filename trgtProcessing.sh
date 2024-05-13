#!/bin/bash

# This shell script downloads and runs TRGT on the file given

# input1: path to input bam
# input2: prefix for naming this sample
# input3: path to reference genome
# input4: path to variant catalog

inputBam=$1
samplePrefix=$2
ref=$3
catalog=$4

wget https://github.com/PacificBiosciences/trgt/releases/download/v0.3.3/trgt-v0.3.3-linux_x86_64.gz
gunzip trgt-v0.3.3-linux_x86_64.gz
mv trgt-v0.3.3-linux_x86_64 trgt
chmod 700 trgt

./trgt --genome ${ref} --repeats ${catalog} --reads ${inputBam} --threads $(nproc) --output-prefix ${samplePrefix}

