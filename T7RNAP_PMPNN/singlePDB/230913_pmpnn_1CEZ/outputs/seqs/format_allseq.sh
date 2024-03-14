#!/bin/bash

# Check if an argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input_filename>"
    exit 1
fi

input_file="$1"
base_name=$(basename "$input_file" .fa) # Extract base name without the .fa extension
formatted_file="${base_name}_formatted.fa"

mkdir -p allfasta

# Remove "T=0.3, " from the input file and save to a formatted file
sed 's/T=0\.3, //g' "$input_file" > "$formatted_file"

# Separate each sequence into a new fasta file under the allfasta folder
awk '{filename = "allfasta/file" int((NR-1)/2); print > filename".txt"}' "$formatted_file"


