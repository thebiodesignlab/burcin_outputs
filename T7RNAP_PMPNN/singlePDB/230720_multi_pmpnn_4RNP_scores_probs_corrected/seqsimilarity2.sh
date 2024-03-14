#!/bin/bash
address="outputs/seqs/pairs/"
# Loop through all the files

find "$address" -type f -name "*.txt" | sort -V | while read -r file; do
  # Extract the samples
  sample1=$(awk -F'=' '/^# 1: / { print $2 }' "$file" | tr -d '_')
  sample2=$(awk -F'=' '/^# 2: sample/ { print $2 }' "$file" | tr -d '_')

  # Check for special case
 if [[ "$sample1" == "4RNP" || -z "$sample1" ]]; then
    sample1="0"
 fi 
  # Extract the similarity score
	similarity=$(awk -F'[()]' '/Similarity:/ { print $2 }' "$file")
  # Print the result
  echo "$sample1 $sample2 $similarity"
  rm $file
done




