#!/bin/bash

conda init bash
source ~/.bashrc
conda activate /camp/home/acarb/.conda/envs/platereader

python3 user_friendly_beta.py
conda deactivate

echo "Done"
