#!/bin/bash
#SBATCH --job-name=epigraph_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=epigraph_test.log

#run your code
echo "Job started at:"
date

echo "Current working directory:"
pwd

echo $SHELL

python3 scripts/epigraph_inference_custom.py --pdb 3CX3 --device cuda --pdb_path test_data/PDB --save_path temp/PDB_Processed --model_path src/EpiGraph/checkpoint --out_path temp/epigraph


echo "Job finished at:"
date