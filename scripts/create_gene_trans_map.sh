#!/bin/bash

transcriptome=$1
output=$2

awk -F " " '/^>/ {print substr($1, 2, index($1, "t")-2) "\t" substr($1, 2)}' $transcriptome > $output 
echo "gene_trans_map were created successfully!"
