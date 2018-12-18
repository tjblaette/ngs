#!/bin/bash

DB="/NGS/known_sites/hg19/biomart_ensembl74ID_to_HUGOsymbol_mapped.txt"

# INPUT files from Limma analysis -> must not contain any "-" in its file name!
IN="$1" #=*_DESeq2results.txt_annotated_woutNA.txt
#ALL OTHER INPUTS = GENE SET COLLECTIONS TO TEST

# Annotate results from DEG analysis with HUGO gene symbols to match gene identifiers in GSEA gene sets
#  $DB contains ENSEMBL IDs in column 1 and HUGO gene symbols in column 2
#  ENSEMBL IDs that could not be mapped to a HUGO gene symbol have already been removed
join -t '	'  <(tail -n +2 $DB | sort -k1,1b) <(tail -n +2 $IN | sort -k1,1b) > ${IN%.txt}_geneSymbols.txt

# Prepare preranked gene list for GSEA
#  Must contain HUGO gene symbols in column 1 and $(sign of logFC) * log(adj.pVal) in column 2
#  Careful: AWK's log function takes the natural logarithm! Use the following for conversion: logb(X) = loga(X)/loga(b) <=> log10(X) = ln(X)/ln(10)
#  Here: $2 = HUGO gene symbol, $5 = logFC, $7 = Wald statistic, $9 = adjusted p-value   
awk -v OFS='\t'	'{if($9 == 0){$9=1.00000e-307} LOG=log($9)/log(10); print($2,LOG)}' ${IN%.txt}_geneSymbols.txt > ${IN%.txt}_forGSEA_combined.rnk
awk -v OFS='\t'	'function sign(v) {return v < 0 ? -1 : 1} {if($9 == 0){$9=1.00000e-307} LOG=log($9)/log(10); print($2,LOG*sign($5))}' ${IN%.txt}_geneSymbols.txt > ${IN%.txt}_forGSEA_directed.rnk
awk -v OFS='\t'	'{print($2,$7)}' ${IN%.txt}_geneSymbols.txt > ${IN%.txt}_forGSEA_stat.rnk
awk -v OFS='\t'	'{print($2,$5)}' ${IN%.txt}_geneSymbols.txt > ${IN%.txt}_forGSEA_lfc.rnk

# Run GSEA once for each gene set
# GSEA will fail if JAVA tmp folder does not exist -> create it now
mkdir -p JAVA_TMP
shift 1
while [[ $# > 0 ]]
do
    set=$1
    java -cp /NGS/gsea/gsea2-2.2.1.jar -Xmx1024m xtools.gsea.GseaPreranked -gmx $set -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ${IN%.txt}_forGSEA_combined.rnk -scoring_scheme weighted -rpt_label "$(basename ${IN})_$(basename ${set})_combined" -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -gui false
    java -cp /NGS/gsea/gsea2-2.2.1.jar -Xmx1024m xtools.gsea.GseaPreranked -gmx $set -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ${IN%.txt}_forGSEA_directed.rnk -scoring_scheme weighted -rpt_label "$(basename ${IN})_$(basename ${set})_directed" -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -gui false
    java -cp /NGS/gsea/gsea2-2.2.1.jar -Xmx1024m xtools.gsea.GseaPreranked -gmx $set -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ${IN%.txt}_forGSEA_stat.rnk -scoring_scheme weighted -rpt_label "$(basename ${IN})_$(basename ${set})_stat" -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -gui false
    java -cp /NGS/gsea/gsea2-2.2.1.jar -Xmx1024m xtools.gsea.GseaPreranked -gmx $set -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk ${IN%.txt}_forGSEA_lfc.rnk -scoring_scheme weighted -rpt_label "$(basename ${IN})_$(basename ${set})_lfc" -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -gui false

# alternative scoring scheme = "classic"    

    shift 1
done



