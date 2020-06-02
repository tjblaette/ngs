#!/bin/bash

####
# T.J.Blätte
# 2015
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
#   REF: Reference genome to extract the flanking sequence from,
#       must match the reference used during alignment and variant
#       calling of the mutations.
#   LEN: Number of flanking bases to annotate on either side of the
#       given mutations.
#
# Output:
#   Input files are annotated in-place and overwritten!
#
####

IN="$1"
REF="$2"
LEN="$3"
TMP=${IN}_tmp


# print date and time of when analysis is started
echo "$(date) running $0"


get_flanking () {
#$1 = chrom
#$2 = pos
#$3 = reference genome fasta
#$4 = number of surrounding bases to print
#$5 = input file to annotate

chr="$1"
pos=$2
ref=$3
length=$4
input="$(basename ${5%.csv})"
#length=$(( 2 * $4 + 1 ))

#make sure pos is not negative - would mean that there is not enough flanking sequence upstream - set to 0 and adjust length to only print what is available
if [ $pos -lt 1 ]
then
	length=$(($length + $pos - 1))
	pos=1
fi

#number of bases per line in reference fasta (50 for GATK's ucsc.hg19 but 60 for the other hg19!)
width=$( tail -2 $ref | wc -L )

#parse non-number-chromosomes (and make sure all others are matched correctly too!)
if [ ! -e ${ref}_chromosomes.txt ]
then
	cut -f1 ${ref}.fai | nl > ${ref}_chromosomes.txt
fi

chr=$(grep -w $chr ${ref}_chromosomes.txt | cut -f1)

#get final chromosome in the reference to check later on that no bases are printed which lie outside the reference
final_chr=$(echo $(tail -n 1 ${ref}_chromosomes.txt | cut -f1) | sed 's/>chr\(.*\)/\1/')

#get the number of lines belonging to each chromosome in the reference fasta
if [ ! -e ${ref}_linesPerChr.txt ]
then
	grep -n '>chr' $ref | cut -f1 -d':' > ${ref}_linesPerChr.txt
	# also append total line number of REF to check later and prevent recalculation
	wc -l $ref | cut -f1 -d' ' >> ${ref}_linesPerChr.txt
fi

#get the region to be printed out
lines_before=$(( $(sed $chr'q;d' ${ref}_linesPerChr.txt) +1 ))
lines_broken=$(( $pos / $width ))
line=$(( $lines_before + $lines_broken ))

start=$(( $pos % $width -1))
#start=$(( $pos % $width - $4))

if [ $start -lt 0 ]
then
	if [ $line -eq $(( $lines_before )) ]
	then
		start=0
	else
		line=$(( $line - 1 ))
		start=$(( $start + $width ))
	fi
fi

end=$(( $start + $length ))


#test if there are enough bases in this row to output full flanking sequences
currentLine=$( sed $line'q;d' $ref )
currentWidth=$( echo $currentLine | wc -L )
if [ $end -gt $currentWidth ]
then
#echo 'too long'

    if [ $line -lt $(tail -n 1 ${ref}_linesPerChr.txt) -a $chr = $final_chr -o $chr != $final_chr -a $line -lt $(( $(sed $(( $chr + 1 ))'q;d' ${ref}_linesPerChr.txt) -1 )) ]
    then
#echo 'have a second line'
	line2=$(( $line + 1 ))
	start2=$(( 1 ))
	end2=$(( $end - $width ))
	end1=$width
#kann ich hier garantieren, dass die nächste Zeile immer reicht??
	echo $currentLine | sed "s/^\(.\{${start}\}\)\(.*\)/\2/g" > ${input}_${1}_${2}_sed_tmp.txt
	sed $line2'q;d' $ref | sed "s/^\(.\{${end2}\}\).*/\1/g" >> ${input}_${1}_${2}_sed_tmp.txt
	
	tr -d '\n' < ${input}_${1}_${2}_sed_tmp.txt | cat
	rm -f ${input}_${1}_${2}_sed_tmp.txt
	echo ''
    else
#echo 'no second line'
        echo $currentLine | sed "s/^\(.\{${start}\}\)\(.*\).*/\2/g"
    fi
else
#echo 'all good'
	#print
	echo $currentLine | sed "s/^\(.\{${start}\}\)\(.\{${length}\}\).*/\2/g"
fi

}

echo "$(head -n1 "$IN"),\"leftFlankingSeq\",\"rightFlankingSeq\"" > $TMP

tail -n +2 $IN | while read line
do
                chrom=$(echo $line | sed 's/"//g' | cut -f1 -d',')
                # indel_length=$(echo $line | sed 's/"//g' | cut -f5 -d',' | grep '-' | sed 's/^-\([ATCG]*\)/\1/' | wc -L)
		indel_begin=$(echo $line | sed 's/"//g' | cut -f2 -d',' )
		indel_end=$(echo $line | sed 's/"//g' | cut -f3 -d',' )
                left=$(( $indel_begin - $LEN ))
                right=$(( $indel_end + 1 ))
		# if it's an insertion (Ref allele == '-'), WT base at variant pos coord must be printed too
            	if [ $(echo "$line" | cut -f4 -d',') == '"-"' ]
            	then
               		paste -d ',' <(echo ${line}) <(echo "\"$(get_flanking $chrom $(( $left + 1 )) $REF $LEN $IN)\",\"$(get_flanking $chrom $right $REF $LEN $IN)\"" | tr 'atcgn' 'ATCGN') >> $TMP
            	else
			paste -d ',' <(echo ${line}) <(echo "\"$(get_flanking $chrom $left $REF $LEN $IN)\",\"$(get_flanking $chrom $right $REF $LEN $IN)\"" | tr 'atcgn' 'ATCGN') >> $TMP
		fi
done

# overwrite input file with tmp output
mv $TMP $IN

