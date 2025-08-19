#!/bin/bash
#SBATCH --job-name=aln_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=aln_test.log

#run your code
echo "Job started at:"
date

echo "Current working directory:"
pwd

echo $SHELL
python3 scripts/standardize_pdbatoms.py --pdb_file test_data/PDB/3CX3.pdb --output_csv temp/3CX3_aln.csv

echo "Job finished at:"
date