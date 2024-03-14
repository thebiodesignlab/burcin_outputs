#!/bin/bash
 
#SBATCH --partition=cpu
#SBATCH --job-name=needle
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --error=err.log
#SBATCH --output=out.log
#SBATCH --mem=1G
#SBATCH --time=01:0:0
#SBATCH --mail-user=brcnacar@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

ml purge > /dev/null 2>&1
module load EMBOSS/6.6.0-foss-2016b
 
# Directory to save the output files
OUTPUT_DIR="pairs"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Directory containing the FASTA files
FASTA_DIR="allfasta"

# List of all FASTA files in the directory
FILES=($(ls ${FASTA_DIR}/*.txt | sort -V))

# Number of files
N=$1
counter=0
for ((i=0; i<N; i++)); do
        for ((j=i+1; j<N; j++)); do
                pair[counter]="$counter $i $j"
                pair_values=(${pair[$counter]})
                k=${pair_values[1]}
                l=${pair_values[2]}
                file1=${FILES[$k]}
                file2=${FILES[$l]}
                let counter++

        done
done
echo $SLURM_ARRAY_TASK_ID
# Get the pair for this array task
ord=$((($2-1)*1000+$SLURM_ARRAY_TASK_ID))
pair_values=(${pair[$ord]})
k=${pair_values[1]}
l=${pair_values[2]}
file1=${FILES[$k]}
file2=${FILES[$l]}

# Run the needle alignment
srun needle -asequence $file1 -bsequence $file2 -gapopen 10.0 -gapextend 0.5 -outfile ${OUTPUT_DIR}/needle_output_${k}_${l}.txt

file=${OUTPUT_DIR}/needle_output_${k}_${l}.txt
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
      echo "$sample1 $sample2 $similarity" > "${OUTPUT_DIR}/extract_${k}_${l}.txt"

  rm $file
