#!/bin/bash


if [ $# -ne 4 ]
then 
    echo "run-tests.sh PROGRAM PARAMS matrices-directory output-file"
    exit
fi

args=("$@")
PROGRAM=${args[0]}
PARAMS=${args[1]}
FILES="${args[2]}/*"
OUTPUT=${args[3]}


for f in $FILES
do
  echo "  <<<<<<< ********************************** >>>>>>>> " | tee -a $OUTPUT
  echo "Processing $f file..." | tee -a $OUTPUT
  
  # take action on each file. $f store current file name
  $PROGRAM $f $PARAMS 2>&1 | tee -a $OUTPUT
  echo "    " | tee -a $OUTPUT
  echo "    " | tee -a $OUTPUT
  echo "    "
  echo "    "
done
