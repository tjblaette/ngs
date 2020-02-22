#!/bin/bash

####
# T.J.Blätte
# 2019
####
#
# Extract some number n of random reads
#       from a given set of paired-end
#       FASTQ files. Save the sampled
#       reads to timestamped files respectively.
#
# Args:
#   IN: Input FASTQ file with forward reads.
#   PAIRED: Paired FASTQ file with reverse reads.
#   N_READS_TO_SAMPLE: Integer value corresponding
#       to the number of reads to be sampled from
#       the two input files.
#
# Output:
#   Out files' prefix depends on that of IN.
#       A suffix is added to denote the number of
#       sampled reads it contains and a timestamp
#       to differentiate distinct subsamples of the
#       same input files. Timestamps include day,
#       hour, minute and second of the subsampling.
#
####

# exit when any error occurs
set -e

IN="$1"
PAIRED="$2"
N_READS_TO_SAMPLE="$3"
OUT_PREFIX="$(echo "$IN" | sed 's/_R.*/_R/')"
OUT_SUFFIX="${N_READS_TO_SAMPLE}_$(date +%d%H%M%S).fastq"


if [ $(( 4 * $N_READS_TO_SAMPLE )) -gt "$(wc -l "$IN" | cut -f1 -d' ')" ]
then
    echo "The input FASTQ files contain fewer reads than you would like to sample!"
else
    READS_TO_SAMPLE="$(sed -n '1~4p' "$IN" | cut -f1 -d' ' | sort -R | head -n $N_READS_TO_SAMPLE)"
    grep -wF -A3 --no-group-separator -f <(echo "$READS_TO_SAMPLE") "$IN" > ${OUT_PREFIX}1_${OUT_SUFFIX}
    grep -wF -A3 --no-group-separator -f <(echo "$READS_TO_SAMPLE") "$PAIRED" > ${OUT_PREFIX}2_${OUT_SUFFIX}
    echo "Successfully sampled ${N_READS_TO_SAMPLE} random reads to ${OUT_PREFIX}?_${OUT_SUFFIX}"
fi
