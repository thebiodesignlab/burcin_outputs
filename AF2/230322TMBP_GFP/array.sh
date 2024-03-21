#!/bin/bash
# This script runs AF2 as a batch job using a GPU and gets 
# folder name from array number given by exec.sh so that
# the AF2 precictions are automated for all fusion sequences.

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
#SBATCH --mail-user=burcin.acar@crick.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

ml purge > /dev/null 2>&1
ml AlphaFold/2.2.2-foss-2021a-CUDA-11.3.1

cd oldAF2
echo $SLURM_ARRAY_TASK_ID
folder="fusion$SLURM_ARRAY_TASK_ID"
echo "$folder"
if [ -d "$folder" ]; then
    for file in "$folder"/*; do
      if [ -f "$file" ]; then
        echo "Reading file: $file"
	srun alphafold --fasta_paths=/camp/lab/debenedictise/home/users/acarb/230322TMBP_GFP/oldAF2/$file --max_template_date=2022-01-01 --run_relax=false --model_preset=monomer --output_dir=/camp/lab/debenedictise/home/users/acarb/230322TMBP_GFP/oldAF2/$folder
      fi
    done
fi
