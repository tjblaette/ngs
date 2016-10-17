#!/bin/bash 

# wrapper script for fusion detection tool SOAPfuse
# takes as input all fastq files to analyze
# creates a list of samples to analyze and creates an output folder for each
# creates symlinks for all input samples in each of these folders
# creates a sample sheet as input to the final SOAPfuse call that analyzes all samples simultaneously


R=1

rm -f soapfuse_sampleSheet.txt
mkdir out

for SAMPLE in "$@"
do
        DIRNAME=$(basename $SAMPLE | cut -d'_' -f1)
        mkdir -p ${DIRNAME}/lib

        ln -s $SAMPLE ${DIRNAME}/lib/run_${R}.fastq

        if [ -z "$DIRLIST" ]
        then
               DIRLIST="$DIRNAME"
               R=2
        else
               if [[ $R == 2 ]]
               then
                      	R=1
                      
                      	READ_LENGTH="$(head "$SAMPLE" | wc -L)"
			echo "$DIRNAME  lib     run     $READ_LENGTH" >> soapfuse_sampleSheet.txt
               else
                      	DIRLIST="$DIRLIST       $DIRNAME"

                      	R=2
               fi
        fi
done

for DIR in $DIRLIST
do
        echo $DIR
done
        
perl /NGS/soapfuse-v1.26/SOAPfuse-RUN.pl -c /NGS/soapfuse-v1.26/config/config.txt -fd . -l soapfuse_sampleSheet.txt -o out

