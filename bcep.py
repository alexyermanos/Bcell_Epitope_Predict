import subprocess
import argparse
import sys
import os
import shutil
import pandas as pd
import glob

def run_bepipred3(name, fasta_dir, temp_dir, pred_model, pdb_dir=None):
    """
    Run Bepipred3 using <fasta_dir>/<name>.fasta.
    If the FASTA does not exist but a PDB is available in pdb_dir, extract it.
    """
    script_path = "scripts/bepipred3_custom.py"

    # Use fasta_dir only if it's a valid string
    if fasta_dir and os.path.isdir(fasta_dir):
        fasta_file = os.path.join(fasta_dir, name + ".fasta")
    else:
        fasta_file = None

    # If FASTA doesn't exist or fasta_dir was not provided, try generating from PDB
    if not fasta_file or not os.path.isfile(fasta_file):
        if pdb_dir is None or not os.path.isdir(pdb_dir):
            print(f"Error: FASTA not found and no valid PDB directory provided.", file=sys.stderr)
            sys.exit(1)

        pdb_file = os.path.join(pdb_dir, name + ".pdb")
        if not os.path.isfile(pdb_file):
            print(f"Error: Neither FASTA ({fasta_file}) nor PDB ({pdb_file}) found.", file=sys.stderr)
            sys.exit(1)

        # Generate temporary FASTA in temp_dir/fasta/
        fasta_tmp_dir = os.path.join(temp_dir, "fasta")
        os.makedirs(fasta_tmp_dir, exist_ok=True)
        fasta_file = os.path.join(fasta_tmp_dir, name + ".fasta")

        try:
            subprocess.run([
                "python3", "-u", "scripts/pdbseqres2fasta.py",
                pdb_file, fasta_file
            ], check=True)
            print(f"Extracted FASTA from {pdb_file} -> {fasta_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error extracting FASTA: {e}", file=sys.stderr)
            sys.exit(1)

    # Now run Bepipred3
    try:
        os.makedirs(os.path.abspath(temp_dir), exist_ok=True)
        subprocess.run([
            "python3", script_path,
            "-i", fasta_file,
            "-o", temp_dir,
            "-pred", pred_model,
            "-add_seq_len"
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Bepipred3: {e}", file=sys.stderr)



def run_discotope3(name, pdb_dir, temp_dir):
    """
    Run Discotope3 using <pdb_dir>/<pdb>.pdb.
    """
    script_path = "src/discotope3_web/discotope3/main.py"
    models_dir = "src/discotope3_web/models"
    pdb_file = os.path.join(pdb_dir, name + ".pdb")

    if not os.path.isfile(pdb_file):
        print(f"Error: PDB file not found at {pdb_file}", file=sys.stderr)
        sys.exit(1)

    try:
        os.makedirs(os.path.abspath(temp_dir), exist_ok=True)
        subprocess.run([
            "python3", script_path,
            "--pdb_or_zip_file", pdb_file,
            "--out_dir", temp_dir,
            "--models_dir", models_dir
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Discotope3: {e}", file=sys.stderr)


def run_epigraph(name, pdb_dir, device, save_path, out_dir):
    """
    Run EpiGraph using <pdb_dir>/<pdb>.
    """
    script_path = "scripts/epigraph_inference_custom.py"
    models_dir = "src/EpiGraph/checkpoint"
    pdb_file = os.path.join(pdb_dir, name + ".pdb")

    if not os.path.isfile(pdb_file):
        print(f"Error: PDB file not found at {pdb_file}", file=sys.stderr)
        sys.exit(1)

    try:
        subprocess.run([
            "python3", script_path,
            "--pdb", name,              # EpiGraph expects PDB ID, not full path
            "--device", device,
            "--pdb_path", pdb_dir,      # directory containing the PDB
            "--save_path", save_path,
            "--model_path", models_dir,
            "--out_path", out_dir,
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running EpiGraph: {e}", file=sys.stderr)
        
    
def process_bepipred3(name, tool_dir):
    """
    Process Bepipred3 output into dataframe keyed on (chain, seqresno, seqresid).
    """
    csv_file = os.path.join(tool_dir, f"bp3_{name}.csv")  # e.g., bp3_3CX3.csv
    if not os.path.isfile(csv_file):
        print(f"Warning: Bepipred3 CSV not found at {csv_file}")
        return None

    df = pd.read_csv(csv_file)

    # Split Accession into structure and chain
    accession_split = df["Accession"].str.split(":", expand=True)
    df["structure"] = accession_split[0]
    df["chain"] = accession_split[1]

    # Residue -> seqresid
    df["seqresid"] = df["Residue"]

    # seqresno: row number within each chain starting from 1
    df["seqresno"] = df.groupby("chain").cumcount() + 1

    # Rename scores
    df = df.rename(columns={
        "BepiPred-3.0 score": "bepipred3_score",
        "BepiPred-3.0 linear epitope score": "bepipred3_epiclass"
    })

    # Keep relevant columns
    df = df[["structure", "chain", "seqresno", "seqresid",
             "bepipred3_score", "bepipred3_epiclass"]]

    return df


def process_discotope3(name, tool_dir):
    """
    Process Discotope3 output into dataframe keyed on (chain, resno, resid).
    Expects tool_dir/<name>/ with chain-specific CSV files:
        e.g., 3CX3_A_discotope3.csv, 3CX3_B_discotope3.csv
    """
    name_dir = os.path.join(tool_dir, name)
    if not os.path.isdir(name_dir):
        print(f"Warning: Discotope3 directory not found at {name_dir}")
        return None

    all_dfs = []

    # Glob all *_discotope3.csv files
    csv_files = glob.glob(os.path.join(name_dir, "*_discotope3.csv"))
    if len(csv_files) == 0:
        print(f"Warning: No Discotope3 CSVs found in {name_dir}")
        return None

    for csv_file in csv_files:
        df = pd.read_csv(csv_file)

        # Extract pdb name from csv filename (before _[chain])
        pdb_base = os.path.basename(csv_file).split("_")[0]

        # Rename columns
        df = df.rename(columns={
            "res_id": "resno",
            "residue": "resid",
            "DiscoTope-3.0_score": "discotope3_score",
            "calibrated_score": "discotope3_calibrated_score",
            "rsa": "discotope3_rsa"
        })

        # Add structure column
        df["structure"] = pdb_base

        # Keep only relevant columns
        df = df[["structure", "chain", "resno", "resid",
                 "discotope3_score", "discotope3_calibrated_score", "discotope3_rsa"]]

        all_dfs.append(df)

    # Concatenate all chains
    discotope3_df = pd.concat(all_dfs, ignore_index=True)
    return discotope3_df


def process_epigraph(pdb, tool_dir):
    """
    Process EpiGraph output into dataframe keyed on (chain, resno, resid).
    """
    csv_file = os.path.join(tool_dir, f"{pdb}.csv")
    if not os.path.isfile(csv_file):
        print(f"Warning: Epigraph file not found at {csv_file}")
        return None 
    
    df = pd.read_csv(csv_file)

    # Split Residue into chain, resid, resno
    residue_split = df["Residue"].str.split(":", expand=True)
    df["chain"] = residue_split[0]
    df["resid"] = residue_split[1]
    df["resno"] = pd.to_numeric(residue_split[2], errors="coerce")

    # Rename columns
    df = df.rename(columns={
        "Score": "epigraph_score",
        "Epitope": "epigraph_epiclass",
        "RSA": "epigraph_rsa"
    })

    # Add required columns
    df["structure"] = df["PDB"]
    df["model"] = 0 #Epigraph by default will only produce outputs for model 0

    # Keep only relevant columns
    df = df[["structure", "model", "chain", "resno", "resid",
             "epigraph_score", "epigraph_epiclass", "epigraph_rsa"]]

    return df
       
        
def standardize_outputs(pdb, pdb_dir, tools, temp_dir, out_dir):
    """
    Standardize outputs and integrate them into a single csv file.
    """
    script_path = "scripts/standardize_pdbatoms.py"
    
    # Step 1: run standardize_pdbatoms
    os.makedirs(out_dir, exist_ok=True)
    pdb_file = os.path.join(pdb_dir, pdb + ".pdb")
    
    if not os.path.isfile(pdb_file):
        print(f"No PDB file found at {pdb_file}. Cannot standardize outputs", file=sys.stderr)
        return  # exit the function early
    
    try:
      subprocess.run([
          "python3", script_path,
          "--pdb_file", pdb_file,
          "--output_dir", os.path.join(temp_dir, "standardized_pdbatoms")
      ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error standardizing outputs: {e}", file=sys.stderr)
    
    # Load standardized pdb atoms as base df
    base_csv = os.path.join(temp_dir, "standardized_pdbatoms", f"{pdb}_standardized.csv")
    if not os.path.isfile(base_csv):
        print(f"Error: standardized file not found at {base_csv}", file=sys.stderr)
        return

    base_df = pd.read_csv(base_csv)
    
    # Step 2: For each tool, load processed outputs
    for tool in tools:
        tool_temp_dir = os.path.join(temp_dir, tool)

        if tool == "bepipred3":
            df_tool = process_bepipred3(pdb, tool_temp_dir)
            base_df = base_df.merge(
                df_tool,
                on=["structure", "chain", "seqresno", "seqresid"],
                how="left"
            )

        elif tool == "discotope3":
            df_tool = process_discotope3(pdb, tool_temp_dir)
            base_df = base_df.merge(
                df_tool,
                on=["structure", "chain", "resno", "resid"],
                how="left"
            )

        elif tool == "epigraph":
            df_tool = process_epigraph(pdb, tool_temp_dir)
            base_df = base_df.merge(
                df_tool,
                on=["structure", "model", "chain", "resno", "resid"],
                how="left"
            )
            
    # Step 3: Save merged df
    merged_csv = os.path.join(out_dir, f"{pdb}_merged.csv")
    base_df.to_csv(merged_csv, index=False)
    print(f"Saved standardized merged output -> {merged_csv}")    
    
def main():
    parser = argparse.ArgumentParser(description="Run B-cell epitope prediction tools.")
    
    parser.add_argument(
        "--tool", 
        choices=["bepipred3", "discotope3", "epigraph"],
        required=True,
        nargs='+',
        help="Select one or more tools to run (bepipred3, discotope3, epigraph)"
    )

    parser.add_argument(
        "--temp_dir", 
        type=str, 
        default="temp",
        help="Temporary directory to store intermediate outputs (default: 'temp')"
    )

    parser.add_argument(
        "--out_dir", 
        type=str, 
        default="output",
        help="Final output directory (default: 'output')"
    )

    parser.add_argument(
        "--pdb_or_fasta",
        type=str,
        required=True,
        help="Base name of the PDB/FASTA file (without extension). Example: '1abc'"
    )

    parser.add_argument(
        "--pdb_dir",
        type=str,
        help="Directory containing PDB files"
    )

    parser.add_argument(
        "--fasta_dir",
        type=str,
        help="Directory containing FASTA files"
    )
    
    parser.add_argument(
        "--epigraph_device", 
        type=str, 
        default="cuda", 
        help="Device to use for EpiGraph (default: 'cuda')"
    )

    parser.add_argument(
        "--bp3pred", 
        type=str, 
        default="vt_pred", 
        choices=['vt_pred', 'mjv_pred'],
        help="Prediction model for Bepipred3 (default: 'vt_pred')"
    )
    
    parser.add_argument(
        "--standardize_outputs",
        type=bool,
        default=True,
        help="Whether to standardize the outputs and integrate them into a single csv (default: True)."
    )

    args = parser.parse_args()

    for tool in args.tool:
        tool_temp_dir = os.path.join(args.temp_dir, tool)
        os.makedirs(tool_temp_dir, exist_ok=True)

        print(f"Running {tool}...")

        if tool == "discotope3":
            run_discotope3(args.pdb_or_fasta, args.pdb_dir, tool_temp_dir)

        elif tool == "epigraph":
            run_epigraph(args.pdb_or_fasta, args.pdb_dir, args.epigraph_device, tool_temp_dir, tool_temp_dir)

        elif tool == "bepipred3":
            run_bepipred3(args.pdb_or_fasta, args.fasta_dir, tool_temp_dir, args.bp3pred, pdb_dir=args.pdb_dir)

    if args.standardize_outputs:
        os.makedirs(args.out_dir, exist_ok=True)
        standardize_outputs(args.pdb_or_fasta, args.pdb_dir, args.tool, args.temp_dir, args.out_dir)

if __name__ == "__main__":
    main()

