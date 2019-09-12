#!/bin/bash

####
# T.J.Blätte
# 2019
####
#
# Takes a reference genome (.fasta) as input and
#       creates several index files needed by one or
#       more pipeline tools. This script must be run
#       once for each reference genome before the
#       pipelines can be run.
#
# Args:
#   REF: Reference genome fasta file.
#
# Output:
#   ${REF}.fai: Reference genome index file.
#   ${REF%.fasta}.dict: Reference sequence dictionary.
#   ${REF}.genomeFile: Genome file, containing in two columns
#       the name and length of each contig.
#   ${REF}.*: several index files created for/ by BWA.
#
####


REF="$1"

REFFAI="${REF}.fai"
REFDICT="${REF%.fasta}.dict"
GENOME_FILE="${REF}.genomeFile"


# create reference index (.fai) for GATK
echo "Preparing reference genome index..."
samtools faidx "$REF"
echo -e "Done!\n"

# create sequence dictionary (.dict) for GATK
echo "Preparing reference genome dictionary for GATK..."
rm "$REFDICT"  # will fail if file exists already
picard CreateSequenceDictionary R="$REF" O="$REFDICT"
echo -e "Done!\n"

# create genome file for BEDTools
echo "Preparing reference genome file for BEDTools..."
cut -f1-2 "$REFFAI" > "${GENOME_FILE}"
echo -e "Done!\n"

# create index files for BWA
echo "Preparing reference genome index files for BWA..."
bwa index -a bwtsw "$REF"
echo -e "Done!\n"
