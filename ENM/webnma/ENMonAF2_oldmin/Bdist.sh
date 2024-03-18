#!/bin/bash
# Simple SLURM sbatch example
#SBATCH --job-name=NMA_oldmin_R0
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:0
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --error=err.log
#SBATCH --output=out.log

ml Anaconda3/2022.05
conda init bash
source ~/.bashrc
conda activate ENM
mkdir -p BDoutput
# Define input files
num_folders=$(find . -mindepth 1 -maxdepth 1 -type d -name "fusion*"| wc -l)
myVar="webnma ca -s -p BDoutput "

for (( x=0; x<$num_folders; x++ )); do
folder="fusion$x/fusion$x/"
num_files=$(find $folder -mindepth 1 -maxdepth 1 -name "*ranked_0.pdb" | wc -l)
echo $num_files
	for (( y=0; y<$num_files; y++ )); do
		file="fusion$x/fusion$x/ranked_$y.pdb"
		file2="fusion$x/fusion$x/f$x""ranked$y.pdb"
		mv $file $file2
  		if [ -f "$file2" ]; then
			myVar+="$file2 "
			echo $file2
      		fi
   	 done
done

echo $myVar
srun $myVar
echo "DONE"
