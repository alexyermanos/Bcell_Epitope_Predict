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
`python bcep.py --pdb 2BIB --pdb_path test_data/PDB --tools discotope3 epigraph bepipred3`

Arguments:
- `--tools` select one or more tools to run from (`discotope3`, `bepipred3`, `epigraph`)
- `--pdb` Local PDB file **without** .pdb extension (required for `discotope3`, `epigraph`) 
- `--pdb_path` Directory where local PDB files are located required for `discotope3`, `epigraph`, alternative input for `bepipred3`) 
- `--fasta` file path to .fasta file (alternative input for `bepipred3`)
- `--out_dir` Final output directory (default: 'output')

Additional arguments:
- `--tmp_dir` folder name to store output of individual tools in their native format (default: `/temp`)
- `--standardize_outputs` This arguments standardizes the outputs of the different tools to generate a single .csv where the scores for each tool are different columns (default is `True`)
- `--bp3pred` Prediction model for Bepipred3 (default: 'vt_pred')
- `--epigraph_device` Device to use for EpiGraph (default: 'cuda', will fallback to `cpu` if unavailable)

### Notes:
- `bepipred3` can take either `--fasta` or `--pdb` as an input:
  - If a `--pdb` is provided, SEQRES records are extracted and used as the sequence input.
  - If running only `bepipred3`, multi-sequence FASTA files are supported.
- If `--standardize_outputs` is set to `False`, no merged .csv is written to `--out_dir`. Instead, raw tool outputs are available in: /[`--temp_dir`]/[`--tool`].