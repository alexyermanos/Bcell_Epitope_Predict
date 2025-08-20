#!/bin/bash
#SBATCH --job-name=dssp_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=dssp_test.log

#run your code
echo "Job started at:"
date

echo "Current working directory:"
pwd

echo $SHELL

python3 dssp.py --pdb "3CX3" --pdb_path "test_data/PDB"

echo "Job finished at:"
date