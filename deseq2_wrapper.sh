#!/bin/bash


# $1 = input file for DESeq2
# $2 = alpha for DEG statistic (default:0.1)
# $3 = DB for gene id conversion (Ensembl to gene symbol)
# $4 = DB for gene id conversion (Ensembl to Entrez gene ID)
# $5 = min absolute LFC tested for (default:0) -> does not take any effect because R function in deseq does not accept variable for lfc parameter

IN=$1
IN_PREFIX="$(basename "${IN%.tsv}")"
DESIGN=${2:-'~ condition'} # design for testing with DESeq2 (must match column names in $1)
ALPHA=${3:-0.1}
LFC=${4:-0.6} # does not take any effect!
BIOMART=${5:-'/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf'}   #was: /NGS/known_sites/human_ensembl_biomart_gene_ID_to_symbol/mart_export_sorted_woutLRG.txt'}
ENTREZ=${6:-'/NGS/known_sites/hg19/biomart_ensembl74ID_to_entrezID_mapped.txt'}

# check if the gene ID to Symbol conversion table exists already
# if it does not, create it now
if [ ! -e "${BIOMART%.gtf}_geneIDtoSymbol.txt" ]
then
	cat "$BIOMART" | cut -f9 | cut -f1,5 -d';' | sed -e 's/[a-z]\+_[a-z]\+//g' -e 's/;/\t/' -e 's/"//g' -e 's/ //g' | grep -v 'level' | sort -k1,1b | uniq > ${BIOMART%.gtf}_geneIDtoSymbol.txt
fi


# RUN DESEQ2 SCRIPT
Rscript $(dirname $0)/deseq2.R $IN "$DESIGN" $ALPHA $LFC


# ANNOTATE GENE COUNTS with gene symbols in addition to ENSEMBL IDs
for COUNTS in "${IN_PREFIX}_DESeq2results_CountsNormalized.txt" "${IN_PREFIX}_DESeq2results_CountsNormalizedTransformed.txt" "${IN_PREFIX}_DESeq2results_CountsNormalizedTransformed_degs.txt"
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
tail -n +2 ${IN_PREFIX}_DESeq2results.txt | sort > ${IN_PREFIX}_DESeq2results_sorted.txt
sed -i 's/"//g' ${IN_PREFIX}_DESeq2results_sorted.txt

echo -n 'geneID	geneSymbol	' > ${IN_PREFIX}_DESeq2results_annotated.txt
head -n 1 ${IN_PREFIX}_DESeq2results.txt | sed -e 's/"//g' >> ${IN_PREFIX}_DESeq2results_annotated.txt
join -t '	' <(sed 's/\.[0-9]*//' "${BIOMART%.gtf}_geneIDtoSymbol.txt") ${IN_PREFIX}_DESeq2results_sorted.txt >> ${IN_PREFIX}_DESeq2results_annotated.txt

#sort -g -k 7,7 ${IN_PREFIX}_DESeq2results_annotated.txt > ${IN_PREFIX}_DESeq2results_sorted.txt
grep -wv 'NA' ${IN_PREFIX}_DESeq2results_annotated.txt >  ${IN_PREFIX}_DESeq2results_annotated_woutNA.txt

# ANNOTATE with entrez gene IDs as well (for SPIA)
join -t '	' <( tail -n +2 ${IN_PREFIX}_DESeq2results_annotated_woutNA.txt) <(tail -n +2 $ENTREZ) > ${IN_PREFIX}_DESeq2results_annotated_woutNA_wEntrezIDs.txt


#rm ${IN_PREFIX}_DESeq2results_annotated.txt
rm -f ${IN_PREFIX}_DESeq2results_sorted.txt
rm -f ${IN_PREFIX}_DESeq2results_CountsNormalized*.txt
#rm -f ${IN_PREFIX}_DESeq2results.txt
