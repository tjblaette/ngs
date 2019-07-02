#!/bin/bash

####
# T.J.Blätte
# 2018
####
#
# Based on a specified FDR cutoff, return the significantly
#       differentially expressed genes shared by one or
#       more DESeq2 analysis result tables.
#
# Args:
#   ALPHA: FDR cutoff for statistical significance of differentially
#       expressed genes.
#   [...]: All other inputs are DESeq2 result tables (*woutNA.txt),
#       from which significantly differentially expressed genes shared
#       by *all* of the provided files are printed to stdout.
#
####


# test that sufficient input was given
# -> at least ALPHA and one DESeq2 output file
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo -e "At least two arguments must be passed:\n\t- ALPHA cutoff to filter DEGs\n\t-INPUT FILE from DESeq2 (*woutNA.txt)"
    exit 1
fi

# test that $1 is numeric and can be used to filter for DEGs with a max FDR
if  [ ! "$1" -eq "$1" ] 2>/dev/null; then
    echo "The first argument must be a number that can be used to filter for DEGs up to a certain max FDR"
    exit 1
fi


# assign input to named variables
ALPHA=$1

# extract DEGs within the first input file passed
# -> replace _ by tab to convert new DESeq2 format to old
#    (new contains ensemblID_geneSymbol in column 1, old contains the same in two columns)
# => necessary to make script compatible with both old and new files and even a mixture of the two
SHARED="$(sed -e 's/_/\t/' -e 's/^\([A-Z0-9]*\)\.[0-9]*/\1/' "$2" | awk -v alpha=$ALPHA '$8 <= alpha' | cut -f1)"

# for each additional input file, select shared DEGs between $SHARED and the additional input file
for i in "${@:3}"
do
    SHARED="$(grep -xf <(echo "$SHARED") <(sed -e 's/_/\t/' -e 's/^\([A-Z0-9]*\)\.[0-9]*/\1/' "$i" | awk -v alpha=$ALPHA '$8 <= alpha' | cut -f1))"
done

# print final shared DEGs
if [ -e "$SHARED" ]
then
    echo "No DEGs present in ALL of the input files passed."
    exit 0
else
    # annotate ensembl IDs with gene symbols
    join -t '  ' <(echo "$SHARED" | sort) <(sed -e 's/_/\t/' -e 's/^\([A-Z0-9]*\)\.[0-9]*/\1/' "$2" | cut -f1-2 | sort)
fi

