#!/bin/bash

####
# T.J.Blätte
# 2015
####
#
# Wrapper to initiate differential gene expression (DGE)
#   analysis with DESeq2.
#
# Args:
#   IN: Design table. A TSV file specifying at least the
#       columns `fileName`, `sampleName` and the condition
#       to test on for differential expression. Each line
#       describes one sample. Filenames must refer to the
#       respective read count files, which are expected to
#       reside within the `intermediate_files` subfolder,
#       within the directory from which the script is run. These
#       count files are those output by HTSeq or STAR,
#       following the respective formatting in the latter case.
#       All of the information provided in this table will be
#       used to annotate the respective output plots.
#   DESIGN: Design formula to test for, something like
#       `~ condition` for an unpaired analysis or
#       `~ patient + condition` for a paired analysis.
#       Provided factors must be defined as columns in $IN.
#   REFERENCE_LEVEL: The reference level defines relative to what
#       fold changes are calculated. Typically, this is something
#       like `untreated`, `control` or similar. The term must
#       be specified in the condition column of the design table
#       that is tested.
#   ALPHA: False discovery rate (FDR) cutoff that defines differentially
#       expressed genes. Defaults to `0.1`.
#   LFC: The minimum log fold change to test for between groups.
#       Defaults to `0.6` (1.5 on a linear / non-log scale). Set it
#       to `0` (1 on a linear / non-log scale) to test for **any**
#       difference between groups.
#   BIOMART: Database file for ensembl gene ID to gene symbol conversion.
#       Defaults to '/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf'
#
# * Optional input files: *
#   replace_sizeFactors.txt: File containing size factors calculated
#       and used by DESeq2 to normalize for different samples' sequencing
#       depth. If a file with this exact name is present, it will be used
#       instead of the size factors calculated by DESeq2 during the analysis
#       of the data in $IN. These files are saved during each analysis,
#       though with a sample-dependent filename. Rename or symlink them to
#       replace size factors in a subsequent analysis.
#   candidates.txt: File containing genes or other features for which a
#       separate set of output files is to be written. IDs must be provided
#       one per line and must match those in the count files of $IN.
#
# Output:
#   Several output files are written, all prefixed with a combination
#   of $IN and $DESIGN:
#       *counts*txt: Tables containing the normalized and
#        log-transformed read counts of all samples.
#       *woutNA.txt: Table containing the test statistics, fold changes,
#        p- and FDR estimates of all genes.
#       *pdf: PDF plots of principal components 1-6 from PCA and heatmaps
#        of different hierarchical clusterings, for differentially expressed
#        genes and increasing subsets of genes with a high coefficient of
#        variation (CV).
#
####


set -e
set -u


IN="$1"
IN_PREFIX="$(basename "${IN%.tsv}")"

DESIGN="$2"
OUT_PREFIX="${IN_PREFIX}_DESeq2$(echo $DESIGN | sed 's/ //g')"

REFERENCE_LEVEL="$3"
ALPHA=${4:-0.1}
LFC=${5:-0.6}

#BIOMART=${6:-'/NGS/known_sites/human_ensembl_biomart_gene_ID_to_symbol/mart_export_sorted_woutLRG.txt'}
BIOMART=${6:-'/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf'}
GENE_ID_TO_SYMBOL_DICT="${BIOMART%.gtf}_geneIDtoSymbol.txt"

# check if the gene ID to Symbol conversion table exists already
# if it does not, create it now
if [ ! -e "$GENE_ID_TO_SYMBOL_DICT" ]
then
    cat "$BIOMART" | cut -f9 | cut -f1,5 -d';' | sed -e 's/[a-z]\+_[a-z]\+//g' -e 's/;/\t/' -e 's/"//g' -e 's/ //g' | grep -v 'level' | sort -k1,1b | uniq > "$GENE_ID_TO_SYMBOL_DICT"
fi


# RUN DESEQ2 SCRIPT
# -> use the script in the same folder as this script
Rscript "$(dirname $0)/deseq2.R" "$IN" "$DESIGN" "$REFERENCE_LEVEL" "$ALPHA" "$LFC" "$OUT_PREFIX" "$GENE_ID_TO_SYMBOL_DICT"
echo ""

# FIX column header of gene-symbol annotated text files
for FILE in "${OUT_PREFIX}_"*txt
do
    sed -i '1s/^\t/geneID_geneSymbol\t/' $FILE
done

# save separate DESeq2 output file with defined p-values only
for FILE in "${OUT_PREFIX}_all.txt"
do
    [ -f "$FILE" ] || continue
    grep -wv 'NA' "$FILE" >  "${FILE%.txt}_woutNA.txt"
done
