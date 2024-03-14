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
  declare -A requeue_count
  while : ; do
      # check if the job is still running
      if [[ -n $last_job_id ]] && squeue -j $last_job_id >/dev/null ; then
	# Extract task IDs of problematic tasks
        problematic_tasks=$(squeue -j $last_job_id | grep "PD" | grep -E "cpu|needle|acarb" | grep "(launch failed requeued held)" | awk '{print $1}' | cut -d_ -f2)

        # If there are problematic tasks, requeue them
        for task_id in $problematic_tasks; do
            # Increment requeue count for this task_id
            requeue_count[$task_id]=$(( ${requeue_count[$task_id]} + 1 ))

            # If requeued more than 10 times, break the loop
            if [ ${requeue_count[$task_id]} -gt 10 ]; then
                echo "Task $last_job_id"_"$task_id has been requeued 10 times and is still not running. Breaking."
                break 2
            fi

	    echo "Requeuing problematic task $last_job_id"_"$task_id."
            scontrol requeue $last_job_id"_"$task_id
        done
	sleep 30
	echo "sleep"
      else
        # if the job is not running, break the loop and start the next job
        echo "break"
	break
      fi
  done
  unset requeue_count
done
ls pairs/extract_*.txt | sort -V | xargs cat >> extract_data.txt

echo "completed"
