#!/bin/bash
#SBATCH --job-name=bcep_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=bcep_test.log

#run your code
echo "Job started at:"
date

echo "Current working directory:"
pwd

echo $SHELL
python3 bcep.py --fasta "test_data/pdb_protein_sequences_part_2.fasta" --pdb "3CX3" --pdb_path "test_data/PDB" --temp_dir "temp" --out_dir "test_output" --tool bepipred3 epigraph discotope3 

echo "Job finished at:"
date