#!/bin/bash

# $0 = full path of this script itself
JAR_DIR=$( dirname "$(readlink -e "$0")" )

# find all JAR files in the current directory --> use directory of where this script is located, not PWD!
JARS=$(find "$JAR_DIR" -maxdepth 1 -type f -name "*.jar")

if [[ $(echo "$JARS" | wc -l) -gt 1 ]]
then
    echo "more than 1 JAR file found!"
else
    #echo "found exactly 1 JAR file!"
    # execute the JAR file and pass all parameters forward
    java -jar "$JARS" "$@"
fi
