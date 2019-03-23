#!/bin/bash

####
# T.J.Blätte
# 2015
####
#
# Collects all of Picard's metric files from
#       a specified folder, containing
#       alignment and duplicate statistics,
#       and saves them to a single TSV file.
#
# Args:
#   DIR: Name of / Path to the folder to
#       collect Picard files from.
#
# Output:
#   $DIR/bamStat_summary.tsv: containing all
#       of the collected and merged statistics.
#
####


DIR=${1:-'.'}
IN="*metrics.txt"
OUT="${DIR}/bamStat_summary.tsv"

# test whether DIR exists, abort if it does not
if [ ! -d "$DIR" ]
then
    echo "Directory \"${DIR}\" not found"
    echo "Aborting!"
    exit
fi


# test whether DIR contains stat files, abort it it does not
if [ -n "$(find "$DIR" -maxdepth 1 -path $IN -type f)" ]
then
    echo "Processing stat files in \"${DIR}\""
    echo "Procesing..."
    echo ""

    # insert header
    echo -e "File name\t# reads unmapped total\t# reads mapped total\t# paired reads mapped\t# duplicate reads total\t# paired duplicate reads\t# paired optical duplicate reads\t#percent duplicate reads total" > "$OUT"

    # extract mapping and duplicate statistics from all metrics.txt-files in $DIR
    for file in ${DIR}/$IN
    do
        echo "$(basename "$file")"
        echo -ne "$(basename "$file")\t" >> "$OUT"
        head  -n 8 "$file" | tail -n 1 | sed 's/[a-zA-Z]//g' | awk -v FS="\t" -v OFS="\t" '{print $5,$2+2*$3,2*$3,$6+2*$7,2*$7,2*$8,$9}' >> "$OUT"
    done
    echo ""
    echo "Done!"
else
    echo "No stat files found in \"${DIR}\""
    echo "Aborting!"
fi

