#!/bin/bash

# given a threshold ALPHA and one or multiple DESeq2 output files (*woutNA.txt), output DEGs (based on ALPHA) shared by all DESeq2 output files given
# --> if only one file is given, output DEGs of that file
# --> if multiple files are given, output DEGs present in all of these

# test that sufficient input was given (at least ALPHA and one DESeq2 output file)
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
SHARED="$(awk -v alpha=$ALPHA '$8 < alpha' "$2" | cut -f1)"

# for each additional input file, select shared DEGs between $SHARED and the additional input file
for i in "${@:3}"
do
    SHARED="$(grep -xf <(echo "$SHARED") <(awk -v alpha=$ALPHA '$8 < alpha' "$i" | cut -f1))"
    #grep -xf <(echo "$SHARED") <(awk -v alpha=$ALPHA '$8 < alpha' "$i" | cut -f1)
done

# print final shared DEGs
if [ -e "$SHARED" ]
then
    echo "No DEGs present in ALL of the input files passed."
    exit 0
else
    # annotate ensembl IDs with gene symbols
    join -t '	' <(echo "$SHARED" | sort) <(cut -f1-2 "$2" | sort)
fi

