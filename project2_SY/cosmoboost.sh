#!/bin/bash
#SBATCH --job-name=cosmoboost
#SBATCH --output=output.txt
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=24GB
#SBATCH --mail-user=maamari@usc.edu
#SBATCH --mail-type=ALL

cd /home/rcf-proj3/ep4/maamari/CosmoBoost
python anlJeong.py
