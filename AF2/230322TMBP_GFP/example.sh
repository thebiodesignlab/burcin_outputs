#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --job-name=TMBPnewAF2
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --error=err.log
#SBATCH --output=out.log
#SBATCH --gres=gpu:1
#SBATCH --mem=0
#SBATCH --time=21-0:0:0

num_folders=$(find oldAF2 -mindepth 1 -maxdepth 1 -type d | wc -l)
#SBATCH --array=1-$num_folders
#SBATCH --mail-user=brcnacar@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

ml purge > /dev/null 2>&1
ml AlphaFold/2.2.2-foss-2021a-CUDA-11.3.1

echo "$num_folders"
cd oldAF2
for (( x=0; x<=$num_folders; x++ )); do
  folder="fusion$x"
  if [ -d "$folder" ]; then
    for file in "$folder"/*; do
      if [ -f "$file" ]; then
        echo "Reading file: $file"
	srun alphafold --fasta_paths=/camp/lab/debenedictise/home/users/acarb/230322TMBP_GFP/oldAF2/$file --max_template_date=2022-01-01 --run_relax=false --model_preset=monomer --output_dir=/camp/lab/debenedictise/home/users/acarb/230322TMBP_GFP/oldAF2/$folder
      fi
    done
  fi
done
