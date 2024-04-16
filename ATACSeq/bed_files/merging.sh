#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <BED file 1> <BED file 2> <Output file name>"
    exit 1
fi

# Assign input arguments to variables
BED_FILE1=$1
BED_FILE2=$2
OUTPUT_FILE=$3

# Concatenate the two BED files
cat ${BED_FILE1} ${BED_FILE2} > combined.bed

# Sort the combined BED file
sortBed -i combined.bed > sorted_combined.bed

# Merge the sorted BED file
bedtools merge -i sorted_combined.bed -c 4 -o mean > ${OUTPUT_FILE}


echo "Merged BED file created: ${OUTPUT_FILE}"
