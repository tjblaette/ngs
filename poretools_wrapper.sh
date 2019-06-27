#!/bin/bash

####
# T.J.Blätte
# 2016
####
#
# Wrapper script for poretools to convert Oxford nanopore FAST5
#       files to FASTQ format. Several other functions are called
#       to generate some statistics on the data.
#
# Args:
#   DIR: Folder containing the FAST5 data to be converted.
#
# Output:
#   Several output files, all prefixed with the input directory's
#       basename, including FASTQ, PDF and various other text files.
#
####


DIR="$1"
DIRNAME="$(basename "$DIR")"

# create fastq files from fast5
poretools fastq $DIR > ${DIRNAME}.fastq

# plot sequencing yield over time
poretools yield_plot --plot-type reads  --saveas ${DIRNAME}_reads.pdf $DIR
poretools yield_plot --plot-type basepairs  --saveas ${DIRNAME}_basepairs.pdf $DIR

# plot signal over time
#poretools squiggle --saveas ${DIRNAME}_signal.pdf $DIR

# report the longest read
poretools winner $DIR > ${DIRNAME}_longest.txt

# report read size statistics
poretools stats $DIR > ${DIRNAME}_readSize.tsv
poretools hist --saveas ${DIRNAME}_readSize.pdf $DIR

# report base composition statistics
poretools nucdist $DIR > ${DIRNAME}_baseComp.tsv

# report BQS statistics
poretools qualdist $DIR > ${DIRNAME}_bqs.tsv
poretools qualpos --saveas ${DIRNAME}_bqs.pdf $DIR

# report pore performance
poretools occupancy --saveas ${DIRNAME}_occupancy.pdf $DIR
