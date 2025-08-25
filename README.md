# BCEP-pipeline

The pipeline takes as input fasta protein sequences and/or pdb files, which are passed into various pretrained models and outputs per-residue epitope propensity scores as a csv.

## Set up environment

- Create a new environment:
`conda create -n BCEP`

- Activate the environment:
`conda activate BCEP`

### Clone different tools into src
```
cd src
git clone https://github.com/UberClifford/BepiPred3.0-Predictor
git clone https://github.com/Magnushhoie/discotope3_web/
git clone https://github.com/sj584/EpiGraph

cd discotope3_web
unzip models.zip
```

## To install necessary libraries:
```
conda install --file conda_requirements.txt -c conda-forge -c bioconda -c pyg -y
pip install -r pip_requirements.txt
```

## Run the pipeline:
The main script is  `bcep.py`.
### Example usage:
`python bcep.py --pdb_or_fasta 2BIB --pdb_dir test_data/PDB --tools discotope3 epigraph bepipred3`

Arguments:
- `--tools` select one or more tools to run from (`discotope3`, `bepipred3`, `epigraph`)
- `--pdb_or_fasta` local PDB/fasta file **without** .pdb/.fasta extension (required) 
- `--pdb_dir` directory where local PDB files are located (required for `discotope3`, `epigraph`, alternative input for `bepipred3`) 
- `--fasta_dir` directory where local FASTA files are located (alternative input for `bepipred3`) 
- `--out_dir` final output directory (default: 'output')

Additional arguments:
- `--tmp_dir` folder name to store output of individual tools in their native format (default: `/temp`)
- `--standardize_outputs` This arguments standardizes the outputs of the different tools to generate a single .csv where the scores for each tool are different columns (default is `True`)
- `--bp3pred` Prediction model for Bepipred3 (default: 'vt_pred')
- `--epigraph_device` Device to use for EpiGraph (default: `cuda`, will fallback to `cpu` if unavailable)

### Notes:
- `bepipred3` can take either `--pdb_dir` or `--fasta_dir` as an input:
  - If only `--pdb_dir` is provided, SEQRES records are extracted and used as the sequence input.
  - For `--standardize_outputs` its important that when using `--fasta_dir` [pdb_or_fasta].fasta, the header inside the .fasta and [pdb_or_fasta].pdb match exactly, otherwise the outputs can't be integrated with the other tools.  
  - If running only `bepipred3`, multi-sequence FASTA files are supported, but requires --standardize_outputs to be set to False
- If `--standardize_outputs` is set to `False`, no merged .csv is written to `--out_dir`. Instead, raw tool outputs are available in: /[`tmp_dir`]/[`tool`].

### Not yet tested/implemented
- Intended use of `--pdb_dir` and `--fasta_dir`is to allow for processing of multiple .pdb or .fasta files at the same time. Currently, the user has to run bcep.py multiple times. 