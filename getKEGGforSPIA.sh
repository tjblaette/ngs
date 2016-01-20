#!/bin/bash
# script to download present KEGG pathway KGML files to create SPIA data object

ORG=$1 #organism, hsa for human, mmu for mouse


mkdir -p /NGS/KEGG/${ORG}
cd /NGS/KEGG/${ORG}


# get list of pathways for $ORG
wget "http://rest.kegg.jp/link/${ORG}/pathway"

cut -f1 pathway | sort | uniq > pathway_uniq

# get all of these pathways in kgml format and give them meaningful names
while read PATHWAY 
do
    OUT_FILE=$(echo $PATHWAY | cut -f2 -d':')
    wget "http://rest.kegg.jp/get/${PATHWAY}/kgml"
    mv kgml ${OUT_FILE}.xml
done < pathway_uniq

echo "library(\"SPIA\")" > getKEGGforSPIA.R
echo "makeSPIAdata(kgml.path=\"/NGS/KEGG/${ORG}\",organism=\"${ORG}\",out.path=\"/NGS/KEGG/${ORG}\")" >> getKEGGforSPIA.R
Rscript getKEGGforSPIA.R


