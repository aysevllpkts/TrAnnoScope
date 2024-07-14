#!/bin/bash
# examine the number of transcripts that were assembled that appear to be full-length or nearly full-length.

# Store the diamond output file name
diamond_output_file=$1

# Initialize variables
total_hits=$(wc -l < "$diamond_output_file")
cumulative_count=0

# Print table header
echo -e "perc_cov_bin\tcount_in_bin\t>cumulative_bin"

# Loop through percentages from 100 to 0 with a step of -10
for (( i=100; i>=0; i-=10 )); do
    # Calculate the lower and upper bounds for the percentage range
    lower_bound=$((i - 10))
    upper_bound=$((i - 1))

    # Count the hits within the percentage range
    count=$(awk -v lower="$lower_bound" -v upper="$upper_bound" '$NF >= lower && $NF <= upper {count++} END {print count}' "$diamond_output_file")

    # Update cumulative count
    cumulative_count=$((cumulative_count + count))

    # Print the results
    echo -e "[$i-$((i-10)))\t$count\t$cumulative_count"
done
