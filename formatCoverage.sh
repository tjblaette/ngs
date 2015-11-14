#!/bin/bash

IN=$1 #input file is output of BEDTools coverageBed function with -hist option
FASTQ1=$2
FASTQ2=$3
# all other inputs are coverage thresholds of interest


TOTAL_BASES_FASTQ1=$(awk -v LINE=1 -v TOTAL_BASES=0 '{if (LINE == 0) {TOTAL_BASES=TOTAL_BASES+length($0); LINE=3;} else  {LINE=LINE-1}} END {print TOTAL_BASES}' $FASTQ1)
TOTAL_BASES_FASTQ2=$(awk -v LINE=1 -v TOTAL_BASES=0 '{if (LINE == 0) {TOTAL_BASES=TOTAL_BASES+length($0); LINE=3;} else  {LINE=LINE-1}} END {print TOTAL_BASES}' $FASTQ2)

head -n 1 $IN | awk -v FASTQ1=$TOTAL_BASES_FASTQ1 -v FASTQ2=$TOTAL_BASES_FASTQ2 '{print "reads on target: "$8" of "FASTQ1+FASTQ2 " initial reads in fasta ("100* $8/(FASTQ1+FASTQ2)"%)"}'
head -n 1 $IN | awk '{print "on average "$8/$4" reads per bp on the target"}'

for CUTOFF in "${@:4}"
do
	# for all cutoffs except 0, cutoff of minimum depth of interest -> for 0, bases with exactly no coverage is of interest
	if [[ $CUTOFF == 0 ]]
	then
		# 0 coverage proportion is extracted from the original BEDTools column -> conversion to % has not yet been done -> multiply with 100
		awk -v OFS='\t' -v CUTOFF=$CUTOFF '$2==0 {print $3" bp had no coverage ("$5*100"%)"} END {print "0 bp had no coverage (0%)"}' $IN | head -n 1
	elif [[ $CUTOFF == 1 ]]
	then	
		# 1x coverage requires different out "at least 1 read" instead of "at least 1 reads"
		awk -v OFS='\t' -v CUTOFF=$CUTOFF '$2 >= CUTOFF {print $6" bp had at least 1 read ("$7"%)"} END {print "0 bp had at least 1 read (0%)"}' $IN | head -n 1
	else
		awk -v OFS='\t' -v CUTOFF=$CUTOFF '$2 >= CUTOFF {print $6" bp had at least "CUTOFF" reads ("$7"%)"} END {print "0 bp had at least "CUTOFF" reads (0%)"}' $IN | head -n 1
	fi
done

# reads on target: 4350038 of 159936615 initial reads in fasta (2.71985%)
#on average 33.2595 reads per bp on the target
