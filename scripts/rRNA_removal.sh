#!/bin/bash
# get transcript IDS for rRNA hits.

blastn=$1
fasta=$2
output_qseqid=$3
#output_fasta=$4


# Check if the pattern is present in the blastn file
if grep -q -E "rRNA,|, rRNA" $blastn; then
    # If the pattern is found, extract the sequences
    cat $blastn | grep -E "rRNA,|, rRNA" | cut -f 1 | sort | uniq > $output_qseqid
    #seqkit grep -n -f $output_qseqid $fasta -o $output_fasta 
else
    # If the pattern is not found:
    echo "No rRNA fragments were found!"
    touch $output_qseqid 
fi