#!/bin/bash
#SBATCH --job-name=dt3_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=dt3_test.log

#run your code
echo "Job started at:"
date

echo "Current working directory:"
pwd

echo $SHELL
python3 src/discotope3_web/discotope3/main.py --pdb_or_zip_file "test_data/2BIB.pdb" --out_dir "test_data/test_output" --models_dir "src/discotope3_web/models"

echo "Job finished at:"
date