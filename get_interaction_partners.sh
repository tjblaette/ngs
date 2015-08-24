#!/bin/bash

#script takes a list of genes and retriebes known and predicted protein-protein interaction partners using the STRING db

#$1 = gene list

GENES=$1
MINSCORE=${2:-900}
MAXPARTNERS=${3:-10}
DBFOLDER=${4:-'/media/data4/NGS/stringdb'}

#delete file to which results will be appended
rm -f ${GENES}_interactionPartners_9606.v10.STRING.txt


#get human interactions only based on NCBI taxonomy ID (9606 for human)
if ! [ -f ${DBFOLDER}/9606.protein.aliases.v10.txt ]
then
	sed -n '/^9606\./p' ${DBFOLDER}/protein.aliases.v10.txt > ${DBFOLDER}/9606.protein.aliases.v10.txt
fi

#restrict aliases to those of BioMart to get unique entries
if ! [ -f ${DBFOLDER}/biomart.9606.protein.aliases.v10.txt ]
then
	grep 'BioMart' ${DBFOLDER}/9606.protein.aliases.v10.txt > ${DBFOLDER}/biomart.9606.protein.aliases.v10.txt
fi

#retrieve entries for input genes in the alias file
#this enables the translation of gene to ensemble protein ID
while read line
do 
gene=$(echo $line |  sed 's/"//g')
	grep -w "$gene" ${DBFOLDER}/biomart.9606.protein.aliases.v10.txt | sed -n "/\t$gene\t/p" | cut -f1 >> ${GENES}_biomart.9606.protein.aliases.v10.txt 
done < $GENES 

#take ensemble IDs of the input genes and search for protein interactions with a certain minimum score and a maximum number of interaction partners to report
#cut to isolate the interaction partners' ensemble protein ID only
while read id
do 
	awk -v id="$id" -v MINSCORE="$MINSCORE" '($1 == id) && ($3 >= MINSCORE) {print $0}' ${DBFOLDER}/9606.protein.links.v10.txt | sort  -r -k 3,3 | head -n "$MAXPARTNERS" | cut -f2 -d' '  >> ${GENES}_interactionPartnersENS_9606.protein.links.v10.txt
#	echo $id
done < "${GENES}_biomart.9606.protein.aliases.v10.txt"


#unique interaction partners
sort ${GENES}_interactionPartnersENS_9606.protein.links.v10.txt | uniq > ${GENES}_interactionPartnersENS_uniq_9606.protein.links.v10.txt

#retranslate ensemble ID to gene name -> for all interaction partners but also the original gene input list
while read protein
do
#	echo "$protein"
	grep -w "$protein" ${DBFOLDER}/biomart.9606.protein.aliases.v10.txt | cut -f2  >> ${GENES}_interactionPartners_9606.v10.STRING.txt 	
done < ${GENES}_interactionPartnersENS_uniq_9606.protein.links.v10.txt

sort -o ${GENES}_interactionPartners_9606.v10.STRING.txt ${GENES}_interactionPartners_9606.v10.STRING.txt

#delete temporary files
rm -f ${GENES}_biomart.9606.protein.aliases.v10.txt 
rm -f ${GENES}_interactionPartnersENS_9606.protein.links.v10.txt
rm -f ${GENES}_interactionPartnersENS_uniq_9606.protein.links.v10.txt


