#!/bin/bash
#SBATCH --job-name=bcep_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/bcep_test.log

#run your code
echo "Job started at:"
date

echo "Current working directory:"
pwd

echo $SHELL
python3 bcep.py --pdb_or_fasta "3CX3" --pdb_dir "test_data/PDB" --temp_dir "temp" --out_dir "test_output" --tool epigraph discotope3 bepipred3 

echo "Job finished at:"
date