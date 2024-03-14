#!/bin/bash

# Directory containing the FASTA files
FASTA_DIR="allfasta"

# List of all FASTA files in the directory
FILES=(${FASTA_DIR}/*.txt)

input_file="$1"
base_name=$(basename "$input_file" .fa) # Extract base name without the .fa extension 
first_4_letters="${base_name:0:4}"


# Number of files
N=${#FILES[@]}

nc=$((N * (N - 1) / 2))

counter=22

bs=1000

mx=10000
#printf '\n%.0s' {256000..499500} >> extract_data.txt
# Placeholder for the id of the last job submitted
last_job_id=""

# Loop over batches

for ((i=counter*bs; i<nc; i+=bs)); do
  ((counter++)) # Calculate end of the batch
  end=$((i+bs-1))
  if ((end>nc)); then
    end=nc
  fi
  echo $i
  echo $end
  last_job_id=$(sbatch --array=0-999 pairwise.sh $N $counter $first_4_letters | cut -f4 -d' ')
  echo $last_job_id
  while : ; do
      # check if the job is still running
      if [[ -n $last_job_id ]] && squeue -j $last_job_id >/dev/null ; then
	sleep 30
	echo "sleep"
      else
        # if the job is not running, break the loop and start the next job
        echo "break"
	break
      fi
  done
done
#ls pairs/extract_*.txt | sort -V | xargs cat >> extract_data.txt
find pairs -name 'extract_*.txt' -print0 | sort -zV | xargs -0 cat >> extract_data.txt

echo "completed"
