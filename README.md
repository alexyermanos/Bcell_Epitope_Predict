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
conda install --file conda_requirements.txt -c conda-forge -c bioconda -pyg
pip install -r pip_requirements.txt
```

## Run the pipeline:
Run the script `scripts/bcep.py` with the following arguments:
- `--tools` select one or more tools to run from (`discotope3`, `bepipred3`, `epigraph`)
- `--pdb` file path to a .pdb file (required for `discotope3`, `epigraph`) 
- `--fasta` file path to .fasta file (required for `bepipred3`)
- `--tmp_dir` folder name to store output of individual tools in their native format (default is `/temp`)

Not yet functional: 
- `--out_dir` folder name to store standardized output of the tools