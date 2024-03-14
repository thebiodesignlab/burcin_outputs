#!/bin/bash
# Directory to save the output files
OUTPUT_DIR="pairs"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Directory containing the FASTA files
FASTA_DIR="allfasta"

# List of all FASTA files in the directory
#FILES=(${FASTA_DIR}/*.txt)
FILES=($(ls ${FASTA_DIR}/*.txt | sort -V))

# Number of files
N=$1
counter=0
for ((i=0; i<N; i++)); do
        for ((j=i+1; j<N; j++)); do
                pair[counter]="$counter $i $j"
		#echo ${pair[$counter]}
		pair_values=(${pair[$counter]})
		k=${pair_values[1]}
		l=${pair_values[2]}	
		file1=${FILES[$k]}
		file2=${FILES[$l]}
		#echo $counter
		#echo $file1
		#echo $file2		
                let counter++
		
        done
done
echo $counter

SLURM_ARRAY_TASK_ID=1700
echo $SLURM_ARRAY_TASK_ID
# Get the pair for this array task
pair_values=(${pair[$SLURM_ARRAY_TASK_ID]})
k=${pair_values[1]}
l=${pair_values[2]}
file1=${FILES[$k]}
file2=${FILES[$l]}

file1=${FILES[$k]}
file2=${FILES[$l]}
echo $file1
echo $file2
