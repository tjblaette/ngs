#!/bin/bash

####
# T.J.Blätte
# Feb 2019
####
#
# When placed in the same directory as
#       a single JAR file, executes that
#       JAR file with all arguments passed
#       forward to that JAR.
#
# Args:
#   ... (optional): All arguments are
#       forwarded as is to the JAR file.
#
####


# get the path to the directory that this script resides in
#   --> this is then also the directory that the JAR
#       to be executed must reside in
# $0 = full path of this script itself
JAR_DIR=$( dirname "$(readlink -e "$0")" )

# find all JAR files in that directory
JARS=$(find "$JAR_DIR" -maxdepth 1 -type f -name "*.jar")

# make sure there is only one JAR in this directory
if [[ $(echo "$JARS" | wc -l) -gt 1 ]]
then
    echo "More than 1 JAR file found in ${JAR_DIR}"
    echo "Aborting!"
else
    #echo "found exactly 1 JAR file!"
    # execute the JAR file and pass all parameters forward
    java -jar "$JARS" "$@"
fi
