#!/bin/bash

# Set output directory
output_dir="nt"
mkdir -p "$output_dir"

# URL to the metadata file
metadata_url="https://ftp.ncbi.nlm.nih.gov/blast/db/nt-nucl-metadata.json"

# Fetch the metadata file and extract file URLs
curl -s "$metadata_url" | grep -o '"ftp://[^"]*"' | sed 's/"//g' | sed 's/ftp:/https:/g' | while read -r file_url; do
    file_name=$(basename "$file_url")
    echo "Downloading $file_name..."
    curl -L -o "$output_dir/$file_name" "${file_url}" &&  \
        
    # Extract the .tar.gz file and remove it after extraction
    echo "Extracting $file_name..."
    tar -xf "$output_dir/$file_name" -C "$output_dir" && rm "$output_dir/$file_name"
done

echo "All files downloaded."

