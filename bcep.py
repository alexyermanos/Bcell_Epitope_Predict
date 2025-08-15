import subprocess
import argparse
import sys
import os
import shutil

def run_bepipred3(fasta_file, temp_dir, pred_model):
    script_path = "scripts/bepipred3_custom.py"

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

def run_discotope3(pdb_file, pdb_path, temp_dir):
    script_path = "src/discotope3_web/discotope3/main.py"
    models_dir = "src/discotope3_web/models"
    
    pdb_file_with_ext = pdb_file + ".pdb" #Discotope3 requires the extension
    file_location = os.path.join(pdb_path, pdb_file_with_ext)
    
    try:
        os.makedirs(os.path.abspath(temp_dir), exist_ok=True)
        
        subprocess.run([
            "python3", script_path,
            "--pdb_or_zip_file", file_location,
            "--out_dir", temp_dir,
            "--models_dir", models_dir
        ], check=True)
        
    except subprocess.CalledProcessError as e:
        print(f"Error running Discotope3: {e}", file=sys.stderr)

def run_epigraph(pdb_file, pdb_path, device, save_path, out_dir):
    script_path = "scripts/epigraph_inference_custom.py"
    models_dir = "src/EpiGraph/checkpoint"
    
    try:
        subprocess.run([
            "python3", script_path,
            "--pdb", pdb_file,
            "--device", device,
            "--pdb_path", pdb_path,
            "--save_path", save_path,
            "--model_path", models_dir,
            "--out_path", out_dir,
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running EpiGraph: {e}", file=sys.stderr)

def standardize_outputs(temp_dir, out_dir):
    # Placeholder for future output standardization logic
    pass

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
        "--pdb", 
        type=str, 
        help="PDB file for Discotope3 & EpiGraph (required for those tools)"
    )
    
    parser.add_argument(
        "--pdb_path", 
        type=str, 
        help="Directory where PDB files are located (required for Discotope3 & EpiGraph)"
    )
    
    parser.add_argument(
        "--epigraph_device", 
        type=str, 
        default="cuda", 
        help="Device to use for EpiGraph (default: 'cuda')"
    )

    parser.add_argument(
        "--fasta", 
        type=str, 
        help="FASTA file for Bepipred3 (required unless PDB is provided)"
    )

    parser.add_argument(
        "--bp3pred", 
        type=str, 
        default="vt_pred", 
        choices=['vt_pred', 'mjv_pred'],
        help="Prediction model for Bepipred3 (default: 'vt_pred')"
    )

    args = parser.parse_args()

    for tool in args.tool:
        if tool == "discotope3":
            if not args.pdb:
                print("Error: --pdb is required for Discotope3.", file=sys.stderr)
                sys.exit(1)
            if not args.pdb_path:
                print("Error: --pdb_path is required for Discotope3.", file=sys.stderr)
                sys.exit(1)
            print(f"Running {tool}...")
            run_discotope3(args.pdb, args.pdb_path, args.temp_dir)
        
        elif tool == "epigraph":
            if not args.pdb:
                print("Error: --pdb is required for EpiGraph.", file=sys.stderr)
                sys.exit(1)
            if not args.pdb_path:
                print("Error: --pdb_path is required for EpiGraph.", file=sys.stderr)
                sys.exit(1)
            print(f"Running {tool}...")
            run_epigraph(args.pdb, args.pdb_path, args.epigraph_device, args.temp_dir, args.temp_dir)
        
        elif tool == "bepipred3":
            if not args.fasta:
                if args.pdb:
                    print("No FASTA provided. Extracting SEQRES from PDB...")
                    pdb_file_path = os.path.abspath(os.path.join(args.pdb_path, args.pdb + ".pdb"))
                    if not os.path.isfile(pdb_file_path):
                        print(f"Error: PDB file not found at {pdb_file_path}", file=sys.stderr)
                        sys.exit(1)

                    # Create temp_dir/fasta subdirectory
                    fasta_dir = os.path.abspath(os.path.join(args.temp_dir, "fasta"))
                    os.makedirs(fasta_dir, exist_ok=True)

                    # Path to output FASTA
                    fasta_out = os.path.join(fasta_dir, args.pdb + ".fasta")

                    # Call pdbseqres2fasta.py
                    try:
                        subprocess.run([
                            "python3", "-u","scripts/pdbseqres2fasta.py",
                            pdb_file_path,
                            fasta_out
                        ], check=True)
                        args.fasta = fasta_out  
                    except subprocess.CalledProcessError as e:
                        print(f"Error running pdbseqres2fasta.py: {e}", file=sys.stderr)
                        sys.exit(1)
                else:
                    print("Error: --fasta or --pdb is required for Bepipred3.", file=sys.stderr)
                    sys.exit(1)

            print(f"Running {tool}...")
            run_bepipred3(args.fasta, args.temp_dir, args.bp3pred)

    print("Standardizing outputs...")
    standardize_outputs(args.temp_dir, args.out_dir)


    print("Standardizing outputs...")
    standardize_outputs(args.temp_dir, args.out_dir)

if __name__ == "__main__":
    main()
