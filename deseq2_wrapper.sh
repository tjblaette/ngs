#!/bin/bash


# $1 = input file for DESeq2
# $2 = alpha for DEG statistic (default:0.1)
# $3 = DB for gene id conversion
# $4 = min absolute LFC tested for (default:0) -> does not take any effect because R function in deseq does not accept variable for lfc parameter

IN=$1
ALPHA=${2:-0.1}
BIOMART=${3:-'/NGS/known_sites/hg19/gencode.v19.chr_patch_hapl_scaff.annotation_UCSCcontigs.gtf'}   #was: /NGS/known_sites/human_ensembl_biomart_gene_ID_to_symbol/mart_export_sorted_woutLRG.txt'}
LFC=${4:-0} # does not take any effect!

# check if the gene ID to Symbol conversion table exists already
# if it does not, create it now
if [ ! -e "${BIOMART%.gtf}_geneIDtoSymbol.txt" ]
then
	tail -n +6 "$BIOMART" | cut -f9 | cut -f1,5 -d';' | sed -e 's/[a-z "_]//g' -e 's/;/\t/' | sort -k1,1b | uniq > ${BIOMART%.gtf}_geneIDtoSymbol.txt
fi



Rscript /NGS/myscripts/deseq2.R $IN $ALPHA $LFC

# annotate normalized gene counts with gene symbols in addition to ENSEMBL IDs
# add header
echo -n 'geneID	geneSymbol	' > ${IN}_DESeq2results_countsTable.txt
head -n 1 ${IN}_DESeq2results_CountsTable.txt | sed -e 's/"//g' >> ${IN}_DESeq2results_countsTable.txt
# add annotated gene counts
join -t '	' "${BIOMART%.gtf}_geneIDtoSymbol.txt" <(sort -k1,1b "${IN}_DESeq2results_CountsTable.txt" | sed -e 's/"//g') >> ${IN}_DESeq2results_countsTable.txt



tail -n +2 ${IN}_DESeq2results.txt | sort > ${IN}_DESeq2results_sorted.txt
sed -i 's/"//g' ${IN}_DESeq2results_sorted.txt

echo -n 'geneID	geneSymbol	' > ${IN}_DESeq2results_annotated.txt
head -n 1 ${IN}_DESeq2results.txt | sed -e 's/"//g' >> ${IN}_DESeq2results_annotated.txt
join -t '	' <(sed 's/\.[0-9]*//' "${BIOMART%.gtf}_geneIDtoSymbol.txt") ${IN}_DESeq2results_sorted.txt >> ${IN}_DESeq2results_annotated.txt

#sort -g -k 7,7 ${IN}_DESeq2results_annotated.txt > ${IN}_DESeq2results_sorted.txt
grep -wv 'NA' ${IN}_DESeq2results_annotated.txt >  ${IN}_DESeq2results_annotated_woutNA.txt

#rm ${IN}_DESeq2results_annotated.txt
rm -f ${IN}_DESeq2results_sorted.txt
rm -f ${IN}_DESeq2results_CountsTable.txt
rm -f ${IN}_DESeq2results.txt
