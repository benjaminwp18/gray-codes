#!/bin/bash
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=36:00:00
#SBATCH --output=log.stdout
#SBATCH --mail-user=bwp18@case.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="gray_code_generator"
./build/main

