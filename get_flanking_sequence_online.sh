#!/bin/bash

####
# T.J.Blätte
# 2019
####
#
# Annotate mutations with flanking reference sequence.
#       Input files are expected in CSV format as
#       output by our mutation calling pipelines.
#       Two columns are added to the input file, with
#       left and right flanking sequences of a given length.
#
#       ** Input files are overwritten! **
#
# Args:
#   IN: File with mutations to annotate.
#   FLANKING_SEQ_LENGTH: Number of flanking bases to annotate
#       on either side of the given mutations. Defaults to 15.
#
# Output:
#   Input files are annotated in-place and overwritten!
#
####

IN=$1
FLANKING_SEQ_LENGTH=${2:-15}


# print date and time of when analysis is started
echo "$(date) running $0"


FLANKING_SEQ=$(
        while read line
        do
            RANGE_chr=$(echo "$line" | cut -f1 -d',' | sed 's/"//g')
            RANGE_left=$(echo "$line" | cut -f2 -d',' | sed 's/"//g')
            RANGE_left=$(( RANGE_left - FLANKING_SEQ_LENGTH ))

            # if it's an insertion (Ref allele == '-'), WT base at variant pos coord must be printed too
            if [ $(echo "$line" | cut -f4 -d',') == '"-"' ]
            then
                RANGE_left=$(( RANGE_left +1 ))
            fi

            RANGE_right=$(echo "$line" | cut -f3 -d',' | sed 's/"//g')
            RANGE_right=$(( RANGE_right + FLANKING_SEQ_LENGTH ))

            FLANKING=$(wget -O - -q http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=${RANGE_chr}:${RANGE_left}:${RANGE_right}  | grep -v '^<')
            FLANKING_left=$(echo $FLANKING | grep --only-matching [ACTGactgNn] | head -n $FLANKING_SEQ_LENGTH | tr -d '\n')
            FLANKING_right=$(echo $FLANKING | rev | grep --only-matching [ACTGactgNn] | head -n $FLANKING_SEQ_LENGTH | tr -d '\n' | rev)
            echo  "\"${FLANKING_left}\",\"${FLANKING_right}\"" | tr 'atcgn' 'ATCGN'
        done < <(tail -n +2 "$IN")
        )


echo "$(head -n1 "$IN"),\"leftFlankingSeq\",\"rightFlankingSeq\"" > "${IN}_tmp"
if [ ! -z "$FLANKING_SEQ" ]
then
    paste -d ',' <(tail -n +2 "$IN") <(echo "$FLANKING_SEQ") >> ${IN}_tmp
fi


# overwrite input file with tmp output
mv ${IN}_tmp $IN
