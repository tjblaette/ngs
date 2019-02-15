#!/bin/bash


# $1 = input file for DESeq2, tsv file specifying at least columms fileName, sampleName and the condition to test on
# $2 = design for testing with DESeq2 (must match column names in $1), default: "~ condition"
# $3 = alpha for DEG statistic (independent filtering is used to maximize the numer of genes whose FDR is at most this "alpha"), default: 0.1
# $4 = min absolute log fold change tested for by DESeq2, default: 0.6
# $5 = DB for gene id conversion (Ensembl to gene symbol), default: '/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf'
# $6 = DB for gene id conversion (Ensembl to Entrez gene ID), default: '/NGS/known_sites/hg19/biomart_ensembl74ID_to_entrezID_mapped.txt'

IN=$1
IN_PREFIX="$(basename "${IN%.tsv}")"

DESIGN=${2:-'~ condition'}
OUT_PREFIX="${IN_PREFIX}_DESeq2$(echo $DESIGN | sed 's/ //g')"

ALPHA=${3:-0.1}
LFC=${4:-0.6}
BIOMART=${5:-'/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf'}   #was: /NGS/known_sites/human_ensembl_biomart_gene_ID_to_symbol/mart_export_sorted_woutLRG.txt'}
ENTREZ=${6:-'/NGS/known_sites/hg19/biomart_ensembl74ID_to_entrezID_mapped.txt'}

# check if the gene ID to Symbol conversion table exists already
# if it does not, create it now
if [ ! -e "${BIOMART%.gtf}_geneIDtoSymbol.txt" ]
then
	cat "$BIOMART" | cut -f9 | cut -f1,5 -d';' | sed -e 's/[a-z]\+_[a-z]\+//g' -e 's/;/\t/' -e 's/"//g' -e 's/ //g' | grep -v 'level' | sort -k1,1b | uniq > ${BIOMART%.gtf}_geneIDtoSymbol.txt
fi


# RUN DESEQ2 SCRIPT
# -> use the script in the same folder as this script
Rscript "$(dirname $0)/deseq2.R" "$IN" "$DESIGN" "$ALPHA" "$LFC" "$OUT_PREFIX"


# ANNOTATE GENE COUNTS with gene symbols in addition to ENSEMBL IDs
for COUNTS in "${OUT_PREFIX}_CountsNormalized.txt" "${OUT_PREFIX}_CountsNormalizedTransformed.txt" "${OUT_PREFIX}_CountsNormalizedTransformed_degs.txt"
do
  # define output file name for annotated file
  NEW_FILE=$(echo $COUNTS | sed 's/CountsNormalized/countsNormalized/')

  # add header
  echo -n 'geneID	geneSymbol	' > $NEW_FILE
  head -n 1 "$COUNTS" | sed -e 's/"//g' >> $NEW_FILE
  
  # add annotated gene counts
  join -t '	' <(sed 's/\.[0-9]*//' "${BIOMART%.gtf}_geneIDtoSymbol.txt") <(sort -k1,1b "$COUNTS" | sed -e 's/"//g' -e 's/\.[0-9]*//') >> $NEW_FILE
done


# ANNOTATE DESEQ2 RESULTS with gene symbols
tail -n +2 ${OUT_PREFIX}.txt | sort > ${OUT_PREFIX}_sorted.txt
sed -i 's/"//g' ${OUT_PREFIX}_sorted.txt

echo -n 'geneID	geneSymbol	' > ${OUT_PREFIX}_annotated.txt
head -n 1 ${OUT_PREFIX}.txt | sed -e 's/"//g' >> ${OUT_PREFIX}_annotated.txt
join -t '	' <(sed 's/\.[0-9]*//' "${BIOMART%.gtf}_geneIDtoSymbol.txt") ${OUT_PREFIX}_sorted.txt >> ${OUT_PREFIX}_annotated.txt

#sort -g -k 7,7 ${OUT_PREFIX}_annotated.txt > ${OUT_PREFIX}_sorted.txt
grep -wv 'NA' ${OUT_PREFIX}_annotated.txt >  ${OUT_PREFIX}_annotated_woutNA.txt

# ANNOTATE with entrez gene IDs as well (for SPIA)
join -t '	' <( tail -n +2 ${OUT_PREFIX}_annotated_woutNA.txt) <(tail -n +2 $ENTREZ) > ${OUT_PREFIX}_annotated_woutNA_wEntrezIDs.txt


#rm ${OUT_PREFIX}_annotated.txt
rm -f ${OUT_PREFIX}_sorted.txt
rm -f ${OUT_PREFIX}_CountsNormalized*.txt
#rm -f ${OUT_PREFIX}.txt
