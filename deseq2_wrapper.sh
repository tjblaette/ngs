#!/bin/bash
set -e
set -u

# $1 = input file for DESeq2, tsv file specifying at least columms fileName, sampleName and the condition to test on
# $2 = design for testing with DESeq2 (must match column names in $1), for example: "~ condition" or "~ patient + condition"
# $3 = reference level for the log fold changes calculated by DESeq2, for example "untreated" 
#       -> by default, DESeq2 would choose the first level in alphabetical order 
#       -> to make sure that data is not by mistake misinterpreted, the reference level is set explicitely and must be passed as a parameter by the user
# $4 = alpha for DEG statistic (independent filtering is used to maximize the numer of genes whose FDR is at most this "alpha"), default: 0.1
# $5 = min absolute log fold change tested for by DESeq2, default: 0.6
# $6 = DB for gene id conversion (Ensembl to gene symbol), default: '/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf'
# $7 = DB for gene id conversion (Ensembl to Entrez gene ID), default: '/NGS/known_sites/hg19/biomart_ensembl74ID_to_entrezID_mapped.txt'

IN="$1"
IN_PREFIX="$(basename "${IN%.tsv}")"

DESIGN="$2" #${2:-'~ condition'}
OUT_PREFIX="${IN_PREFIX}_DESeq2$(echo $DESIGN | sed 's/ //g')"

REFERENCE_LEVEL="$3"

ALPHA=${4:-0.1}
LFC=${5:-0.6}
BIOMART=${6:-'/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf'}   #was: /NGS/known_sites/human_ensembl_biomart_gene_ID_to_symbol/mart_export_sorted_woutLRG.txt'}
ENTREZ=${7:-'/NGS/known_sites/hg19/biomart_ensembl74ID_to_entrezID_mapped.txt'}

# check if the gene ID to Symbol conversion table exists already
# if it does not, create it now
if [ ! -e "${BIOMART%.gtf}_geneIDtoSymbol.txt" ]
then
    cat "$BIOMART" | cut -f9 | cut -f1,5 -d';' | sed -e 's/[a-z]\+_[a-z]\+//g' -e 's/;/\t/' -e 's/"//g' -e 's/ //g' | grep -v 'level' | sort -k1,1b | uniq > ${BIOMART%.gtf}_geneIDtoSymbol.txt
fi


# RUN DESEQ2 SCRIPT
# -> use the script in the same folder as this script
Rscript "$(dirname $0)/deseq2.R" "$IN" "$DESIGN" "$REFERENCE_LEVEL" "$ALPHA" "$LFC" "$OUT_PREFIX"
echo ""


# ANNOTATE GENE COUNTS with gene symbols in addition to ENSEMBL IDs
for FILE in "${OUT_PREFIX}.txt" "${OUT_PREFIX}_countsNormalized.txt" "${OUT_PREFIX}_countsNormalizedTransformed.txt" "${OUT_PREFIX}_countsNormalizedTransformed_degs.txt"
do
    # continue only if the file exists
    # --> abort, if there are no DEGs to process
    ([ -f "$FILE" ] && [ "$(grep -c '^ENS' "$FILE")" -gt 0 ]) || continue
    echo "Annotating gene symbols for $FILE"

    # define output file name for annotated file
    ANNOTATED=${FILE%.txt}_annotated.txt

    # add header
    echo -ne 'geneID\tgeneSymbol\t' > $ANNOTATED
    head -n 1 "$FILE" | sed -e 's/"//g' >> $ANNOTATED

    # add and annotated the file
    join -t '	' <(sed 's/\.[0-9]*//' "${BIOMART%.gtf}_geneIDtoSymbol.txt") <(tail -n +2 "$FILE" | sort -k1,1b | sed -e 's/"//g' -e 's/\.[0-9]*//') >> $ANNOTATED

    if [ "$(wc -l "$FILE" | cut -f1 -d' ')" -eq "$(wc -l "$ANNOTATED" | cut -f1 -d' ')" ]
    then
        mv "$ANNOTATED" "$FILE"
    else
        echo "not all features could be annotated - keeping original file"
    fi
done


# save separate DESeq2 output file with defined p-values only
for FILE in "${OUT_PREFIX}.txt" "${OUT_PREFIX}_annotated.txt"
do
    [ -f "$FILE" ] || continue
    grep -wv 'NA' "$FILE" >  "${FILE%.txt}_woutNA.txt"
done


# ANNOTATE with entrez gene IDs as well (for SPIA)
#if [ -f "${OUT_PREFIX}_annotated_woutNA.txt" ]
#then
#    join -t '	' <( tail -n +2 ${OUT_PREFIX}_annotated_woutNA.txt) <(tail -n +2 $ENTREZ) > ${OUT_PREFIX}_annotated_woutNA_wEntrezIDs.txt
#fi
