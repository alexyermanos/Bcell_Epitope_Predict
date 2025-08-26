#!/bin/bash
#SBATCH --job-name=bcep_test
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --array=1-64  # Automatically matches the number of PDBs below

# --- List of PDB/fasta names ---
PDB_LIST=("3CX3" "8QLK" "8QLM" "8QLV" "8QLH" "8QLC" "8QLG" "8A42" "8QM0" "4D0Y" "4HQS" "4HQZ" "2YP6" "6JF1" "5CYB"
          "4H59" "5JJ5" "4HMO" "4HMQ" "4HMP" "6YAB" "6YA3" "6YA4" "6YAG" "6Y9U" "5TVL" "3ZK7" "7L6Z" "4EQ9" "2MVB"
          "5SUO" "5SWA" "5SWB" "4H1X" "2XD2" "2XD3" "2BIB" "1WRA" "4CNL" "4X36" "4IVV" "4Q2W" "7PJ3" "7PJ4" "2WWC"
          "2WW5" "2WWD" "4CP6" "2PMS" "1W9R" "2M6U" "3ZPP" "4MR0" "7F7Y" "4S3L" "2WW8" "2L4O" "2Y1V" "2X9Y" "4QQA"
          "4QQQ" "4E8D" "4E8C" "3E0M" "1W6T" "4CGK" "3ZFJ" "3LFT" "4YZ3")

# Get the PDB/fasta name for this array task
PDB_NAME=${PDB_LIST[$SLURM_ARRAY_TASK_ID-1]}

# Create log files
LOG_FILE="logs/bcep_${PDB_NAME}.log"
exec > "$LOG_FILE" 2>&1

echo "Job started at:"
date

echo "Running for: $PDB_NAME"

python3 bcep.py \
    --pdb_or_fasta "$PDB_NAME" \
    --pdb_dir "test_data/PDB" \
    --temp_dir "temp" \
    --out_dir "test_output" \
    --tool discotope3 epigraph bepipred3 

echo "Job finished at:"
date
