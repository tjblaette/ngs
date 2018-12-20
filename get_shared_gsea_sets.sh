#!/bin/bash

# given a threshold ALPHA and one or multiple GSEA output files (gsea_report*xls), output enriched / depleted gene sets (based on ALPHA) shared by all GSEA output files given
# --> if only one file is given, output gene sets of that file
# --> if multiple files are given, output signficantly enriched / depleted gene sets present in all of these

# test that sufficient input was given (at least ALPHA and one GSEA output file)
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo -e "At least two arguments must be passed:\n\t- ALPHA cutoff to filter gene sets\n\t-INPUT FILE from GSEA (gsea_report*.xls)"
    exit 1
fi

# test that $1 is numeric and can be used to filter for gene sets with a max FDR
if  [ ! "$1" -eq "$1" ] 2>/dev/null; then
    echo "The first argument must be a number that can be used to filter for gene sets up to a certain max FDR"
    exit 1
fi


# assign input to named variables
ALPHA=$1

# extract gene sets within the first input file passed
SHARED="$(awk -v alpha=$ALPHA -v FS='\t' -v OFS='\t' '$8 < alpha' "$2" | cut -f1)"


# for each additional input file, select shared gene sets between $SHARED and the additional input file
for i in "${@:3}"
do
    SHARED="$(grep -xf <(echo "$SHARED") <(awk -v alpha=$ALPHA -v FS='\t' -v OFS='\t' '$8 < alpha' "$i" | cut -f1))"
    #grep -xf <(echo "$SHARED") <(awk -v alpha=$ALPHA -v FS='\t' -v OFS='\t' '$8 < alpha' "$i" | cut -f1)
done

# print final shared gene sets
if [ -e "$SHARED" ]
then
    echo "No gene sets were present in ALL of the input files given."
    exit 0
else
    if [ ! -n "$3" ]
    then
        awk -v FS='\t' -v OFS='\t' -v alpha=$ALPHA 'NR == 1 || $8 < alpha' "$2" | cut -f1,4-10 | sed 's/ /_/g' | column  -t
        exit 0
    else
        echo "$SHARED"
        exit 0
    fi
fi
