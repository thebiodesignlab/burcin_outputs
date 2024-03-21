#!/bin/bash
# This code gives order to start an array of AF2 calculations by running array.sh 
# multiple times (number of fusionX folders - number of different fusion sequences).
# Here the fusionX directories are located under oldAF2.
num_folders=$(find oldAF2 -mindepth 1 -maxdepth 1 -type d | wc -l)
sbatch --array=0-$num_folders array.sh 
