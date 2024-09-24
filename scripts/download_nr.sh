#!/bin/bash

# Set output directory
output_dir="nr"
mkdir -p "$output_dir"

# URL to the metadata file
nr_url="https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"

# Fetch the nr file 
file_name=$(basename "$nr_url")
echo "Downloading $file_name..."
curl -L -o "$output_dir/$file_name" "${nr_url}" 
        
# Extract the .tar.gz file and remove it after extraction
echo "Extracting $file_name..."
tar -xf "$output_dir/$file_name" -C "$output_dir" && rm "$output_dir/$file_name"


echo "NR db were downloaded."
