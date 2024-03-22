# This code extracts the similarity scores.
#!/bin/bash
address="/Volumes/lab-debenedictise/home/users/acarb/T7RNAP_PMPNN/230720_multi_pmpnn_4RNP_scores_probs_corrected/outputs/seqs/pairs/*"

# Loop through all the files
for file in $address; do
  # Extract the samples
  sample1=$(awk -F'=' '/^# 1: / { print $2 }' "$file" | tr -d '_')
  sample2=$(awk -F'=' '/^# 2: sample/ { print $2 }' "$file" | tr -d '_')

  # Check for special case
  if [ "$sample1" == "4RNP" ]; then
    sample1=0
  fi
  
  # Extract the similarity score
  similarity=$(awk '/Similarity:/ { print $3 }' "$file")

  # Print the result
  echo "$sample1 $sample2 $similarity"
done



