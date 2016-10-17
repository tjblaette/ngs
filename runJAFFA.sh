#!/bin/bash

# wrapper script for fusion discovery tool JAFFA
# takes as input gzipped fastq files to analyze
# creates a list of samples
# for each sample it creates an output folder with symlinks to the sample input files
# for-loop then traverses the sample list, enters each folder and analyzes all contained *fastq.gz files


R1="yes"

for SAMPLE in "$@"
do
	DIRNAME=$(basename $SAMPLE | cut -d'_' -f1)
	mkdir -p $DIRNAME
	
	ln -s $SAMPLE $DIRNAME/$(basename $SAMPLE)


	if [ -z "$DIRLIST" ]
	then
		DIRLIST="$DIRNAME"
		R1=""
	else
		if [ -z $R1 ]
		then
			R1="yes"
		else
			DIRLIST="$DIRLIST	$DIRNAME"
			R1=""
		fi
	fi
done

for DIR in $DIRLIST
do
	echo $DIR
	cd $DIR
	/NGS/jaffa-1.06/tools/bin/bpipe run /NGS/jaffa-1.06/JAFFA_direct.groovy *fastq.gz
	cd ..
done
