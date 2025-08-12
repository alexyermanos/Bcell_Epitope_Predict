#!/bin/bash
#SBATCH --job-name=bp3_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=bp3_test.log

#run your code
echo "Job started at:"
date

echo "Current working directory:"
pwd

echo $SHELL
python3 scripts/bepipred3_custom.py -i "test_data/pdb_protein_sequences_part_2.fasta" -o "temp/test_output" -pred "vt_pred" 

echo "Job finished at:"
date