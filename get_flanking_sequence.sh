#!/bin/sh

#$1=file to annotate
#$2=reference genome to use (use the same one as for the variant calling!)
#$3=number of surrounding bases to print

echo "$(date) running $0"

get_flanking () {
#$1 = chrom
#$2 = Pos
#$3 = Reference genome fasta
#$4 = Number of surrounding bases to print
##NOTE: The last base of the left flanking sequence should always equal the reference base at that position --> control for correctness flanking sequences!)

chr="$1"
pos=$2
ref=$3
length=$(($4+1))
#length=$(( 2 * $4 + 1 ))

#make sure pos is not negative - would mean that there is not enough flanking sequence upstream - set to 0 and adjust length to only print what is available
if [ $pos -lt 1 ]
then
	length=$(($length + $pos - 1))
	pos=1
fi

#number of bases per line in reference fasta (50 for GATK's ucsc.hg19 but 60 for the other hg19!)
width=$( head -2 $ref | wc -L ) 

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
	end2=$(( ( $width - $end ) * (-1)))
	end1=$width
#kann ich hier garantieren, dass die nächste Zeile immer reicht??
	echo $currentLine | sed "s/^\(.\{${start}\}\)\(.*\)/\2/g" > ${1}_${2}_sed_tmp.txt
	sed $line2'q;d' $ref | sed "s/^\(.\{${end2}\}\).*/\1/g" >> ${1}_${2}_sed_tmp.txt
	
	tr -d '\n' < ${1}_${2}_sed_tmp.txt | cat
	rm -f ${1}_${2}_sed_tmp.txt
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

head -n 1 $1 > ${1}_tmp
sed -i 's/?$/,left_flanking_seq,right_flanking_seq/' ${1}_tmp

tail -n +2 $1 | while read line
do
                chrom=$(echo $line | cut -f1 -d',') 
                indel_length=$(echo $line | cut -f5 -d',' | grep '-' | sed 's/^-\([ATCG]*\)/\1/' | wc -L)
		indel_begin=$(echo $line | cut -f2 -d',' )
                left=$(( $indel_begin - $3 )) 
                right=$(( $indel_begin +  $indel_length )) 
                echo "${line},$(get_flanking $chrom $left $2 $3),$(get_flanking $chrom $right $2 $3)" >> ${1}_tmp
done

mv ${1}_tmp ${1}_edited

