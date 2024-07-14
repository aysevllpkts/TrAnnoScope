#!/bin/bash
# get transcript IDS for mitochondrial hits of the organism of interest.

blastn=$1
organism=$2
fasta=$3
output_qseqid=$4
#output_fasta=$5

# Check if the pattern is present in the blastn file
if grep -q -E "; mitochondrial|mitochondrion" $blastn; then
    # If the pattern is found, extract the sequences
    cat $blastn | grep -E "; mitochondrial|mitochondrion" | cut -f 1 | sort | uniq > $output_qseqid
    #seqkit grep -n -f $output_qseqid $fasta -o $output_fasta 
else
    # If the pattern is not found:
    echo "No MT fragments were found!"
    touch $output_qseqid
fi