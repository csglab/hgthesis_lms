#!/bin/bash

file1="sim_samplesRM .txt"
file2="run_parameters.txt"
output_file="parameters.txt"

# Check if input files exist
if [ ! -e "$file1" ] || [ ! -e "$file2" ]; then
    echo "Error: Input files do not exist."
    exit 1
fi

# Create the output file with all combinations
while IFS= read -r line1; do
    while IFS= read -r line2; do
        echo "$line1 $line2" >> "$output_file"
    done < "$file2"
done < "$file1"

echo "All combinations have been written to $output_file."

