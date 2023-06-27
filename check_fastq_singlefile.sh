#!/bin/bash

filename="$1"
outputfile="validationresults.txt"

{
## add a line as a file separator between each file output
echo "_________________________"

## add the file name to the output file
echo "File name: $filename"

## followed by the output of fastq_info, which is checking the format of the files

fastq_info -r "$filename" 2>&1

## add an empty line below
echo ""

}  >> "$outputfile"

