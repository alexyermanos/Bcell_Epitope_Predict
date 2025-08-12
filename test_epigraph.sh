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
#python3 src/EpiGraph/inferencev2.py --pdb 1cfi --device cuda --pdb_path test_data/PDB --save_path test_data/PDB_Processed --out_path test_data/test_output

python3 src/epigraph_inference.py --pdb 3CX3 --device cuda --pdb_path test_data/PDB --save_path test_data/PDB_Processed --model_path src/EpiGraph/checkpoint --out_path test_data/test_output

#python3 src/EpiGraph/inference_customPDB.py --pdb "2BIB" --pdb_path "test_data/" --save_path "test_data/PDB_Processed" --model_path src/EpiGraph/checkpoint --out_path "test_data/test_output"

echo "Job finished at:"
date