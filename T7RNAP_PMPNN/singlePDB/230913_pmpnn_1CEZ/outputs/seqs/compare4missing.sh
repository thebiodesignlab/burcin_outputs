#!/bin/bash

# Step 1: Create the full list of expected combinations
for i in {0..999}; do 
  for j in $(seq $((i+1)) 1000); do
    
	echo "$i $j"; 
  done; 
done > full_list.txt

# Step 2: Extract i, j pairs from the file and sort them
awk '{print $1 " " $2}' extract_data.txt | sort -V > existing_pairs.txt

# Step 3: Compare the two lists to find missing pairs
comm -23 full_list.txt existing_pairs.txt

