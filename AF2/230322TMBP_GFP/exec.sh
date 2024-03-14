#!/bin/bash
num_folders=$(find oldAF2 -mindepth 1 -maxdepth 1 -type d | wc -l)
sbatch --array=0-$num_folders array.sh 
