#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --job-name=oldAF2
#SBATCH --ntasks=80
#SBATCH --nodes=1
#SBATCH --error=err.log
#SBATCH --output=out.log
#SBATCH --gres=gpu:4
#SBATCH --mem=0
#SBATCH --time=21-0:0:0
#SBATCH --mail-user=brcnacar@gmail.com
#SBATCH --mail-type=ALL

ml purge > /dev/null 2>&1
ml AlphaFold/2.1.1-fosscuda-2020b
#ml matlab

#matlab -nojvm -nodisplay -r MBP_GFP_DIP > matlab.log

num_folders=$(find oldAF2 -mindepth 1 -maxdepth 1 -type d | wc -l)
echo "$num_folders"
cd oldAF2
for (( x=0; x<=$num_folders; x++ )); do
  folder="fusion$x"
  if [ -d "$folder" ]; then
    for file in "$folder"/*; do
      if [ -f "$file" ]; then
        echo "Reading file: $file"
	alphafold --fasta_paths=/camp/lab/debenedictise/home/users/acarb/230316_MBP_EGFP_GPU/oldAF2/$file --max_template_date=2022-01-01 --model_preset=monomer --output_dir=/camp/lab/debenedictise/home/users/acarb/230316_MBP_EGFP_GPU/oldAF2/$folder
      fi
    done
  fi
done
