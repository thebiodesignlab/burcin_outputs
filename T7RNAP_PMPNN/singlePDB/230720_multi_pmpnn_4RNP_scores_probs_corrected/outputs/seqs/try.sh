#!/bin/bash
# Directory containing the FASTA files
FASTA_DIR="allfasta"

# List of all FASTA files in the directory
FILES=(${FASTA_DIR}/*.txt)

# Number of files
N=$1

# Function to generate unique pairs of files
get_pair() {
    index=$1
    counter=0
    for ((i=0; i<N-1; i++)); do
        for ((j=i+1; j<N; j++)); do
            if ((counter==index)); then
                echo ${FILES[$i]} ${FILES[$j]}
                return
            fi
            let counter++
        done
    done
}
echo $(get_pair 4300)
