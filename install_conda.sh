#!/bin/bash
#SBATCH --job-name=install_conda
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=installation.log

echo "Job started at:"
date

source /hpc/dla_lti/jvanmeenen/miniconda3/etc/profile.d/conda.sh

conda activate TestBCEP

echo "Installing dependencies..."
mamba install --file conda_requirements.txt -c conda-forge -c bioconda -c pyg -y

echo "Job finished at:"
date