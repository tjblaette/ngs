#!/bin/bash

####
# T.J.Blätte
# 2018
####
#
# Based on a specified FDR cutoff, returns the significant gene sets
#       shared by one or more GSEA result tables.
#
# Args:
#   ALPHA: FDR cutoff for statistical significance of gene set
#       enrichment.
#   RANK: Ranking metric used for the GSEA preranked input gene lists.
#       Though previous versions supported various rankings, currently
#       gsea.sh will generate only those based on the log fold change.
#       For these analyses, pass 'lfc'.
#   [...]: All other inputs are GSEA result table *prefixes*. These
#       prefixes must make
#
#       ${PREFIX}*${GENE_SET_COLLECTION}*${RANKING_METRIC}*/gsea_report_*${CHANGE}*.xls
#
#       unique - they are used to collect all the GSEA result tables
#       from which shared significant gene sets are gathered.
#       For each combination of gene-set-collection, ranking-metric and
#       direction-of-enrichment (enrichment vs depletion), the
#       significant gene sets shared by files of *all* of the provided
#       file prefixes are returned and printed to stdout.
####


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
