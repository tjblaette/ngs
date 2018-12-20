#!/bin/bash

# given a threshold ALPHA, the ranking metric RANK used for GSEA preranked input gene lists and one or multiple GSEA output file prefixes (PATH/TO/GSEA/OUTPUT/FOLDER/DESEQ2-OUTPUT_), output enriched / depleted gene sets (based on ALPHA) shared by all GSEA output files with the provided prefix that were ranked by RANK
# --> if only one file is given, output gene sets of that file
# --> if multiple files are given, output signficantly enriched / depleted gene sets present in all of these

# test that sufficient input was given (at least ALPHA, RANK and one GSEA output file prefix)
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo -e "At least three arguments must be passed:\n\t- ALPHA cutoff to filter gene sets\n\t-RANKING METRIC used for GSEA analysis\n\t-INPUT FILE from GSEA (gsea_report*.xls)"
    exit 1
fi

# test that $1 is numeric and can be used to filter for gene sets with a max FDR
if  [ ! "$1" -eq "$1" ] 2>/dev/null; then
    echo "The first argument must be a number that can be used to filter for gene sets up to a certain max FDR"
    exit 1
fi

# assign input to named variables
ALPHA=$1  # 0.25
RANK=$2  # lfc

# loop over combinations of MSIGDB gene sets tested by GSEA and relative enrichment and, for all files that exist with respective details, determine and output overlapping / shared gene sets with max FDR of ALPHA
for SET in "c1" "c2" "c6" "h"
do
    for CHANGE in "pos" "neg"
    do
        echo "${SET} - ${CHANGE}"
        echo "-------------------"
        SAMPLES=""
        for i in ${@:3}
        do
            for file in ${i}*${SET}*${RANK}*/gsea_report_*${CHANGE}*.xls
            do
                if [ -e "$file" ]
                then
                    SAMPLES="$(echo "${SAMPLES} ${file}")"
                fi
            done
        done
            
        if [ -n "$SAMPLES" ]
        then
            get_shared_gsea_sets.sh $ALPHA $SAMPLES
        else
            echo "No samples to process"
        fi
        echo "==================="
    done
done
