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

# download TRGT
wget https://github.com/PacificBiosciences/trgt/releases/download/v0.3.3/trgt-v0.3.3-linux_x86_64.gz
gunzip trgt-v0.3.3-linux_x86_64.gz
mv trgt-v0.3.3-linux_x86_64 trgt
chmod 700 trgt

# prepare to run segmentation script
pip install pandas Levenshtein regex

# genotype tandem repeats with TRGT
./trgt --genome ${ref} --repeats ${catalog} --reads ${inputBam} --threads $(nproc) --output-prefix ${samplePrefix}

# parse TRGT output specifically for the FGF14 locus
python segmentTRGTAlleles.py --sample=${samplePrefix}

