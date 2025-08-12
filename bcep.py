import subprocess
import argparse
import sys
import os
import shutil

# Function to run Bepipred3
def run_bepipred3(fasta_file, temp_dir, pred_model):
    script_path = "src/BepiPred3.0-Predictor/bepipred3_CLI.py"

    try:
        # Ensure the temp directory exists
        os.makedirs(temp_dir, exist_ok=True)
        
        # Run Bepipred3
        subprocess.run([
            "python3", script_path,
            "-i", fasta_file,
            "-o", temp_dir,  # Output to temp directory
            "-pred", pred_model
        ], check=True)
        
    except subprocess.CalledProcessError as e:
        print(f"Error running Bepipred3: {e}")

# Function to run Discotope3
def run_discotope3(pdb_file, pdb_path, temp_dir):
    script_path = "src/discotope3_web/discotope3/main.py"
    models_dir = "src/discotope3_web/models"
    
    file_location = os.path.join(pdb_path, pdb_file)
    
    try:
        # Ensure the temp directory exists
        os.makedirs(temp_dir, exist_ok=True)
        
        # Run Discotope3
        subprocess.run([
            "python3", script_path,
            "--pdb_or_zip_file", file_location,
            "--out_dir", temp_dir,  # Output to temp directory
            "--models_dir", models_dir
        ], check=True)
        
    except subprocess.CalledProcessError as e:
        print(f"Error running Discotope3: {e}")

# Function to run EpiGraph        
def run_epigraph(pdb_file, pdb_path, device, save_path, out_dir):
    script_path = "src/epigraph_inference.py"
    models_dir = "src/EpiGraph/checkpoint"
    
    try:
        subprocess.run([
            "python3", script_path,
            "--pdb", pdb_file,
            "--device", device,
            "--pdb_path", pdb_path,
            "--save_path", save_path,
            "--model_path", models_dir,
            "--out_path", out_dir
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running EpiGraph: {e}")

# FIXME: Function to standardize outputs depending on tool (left empty for now)
def standardize_outputs(temp_dir, out_dir):
    # Placeholder for future output standardization logic
    pass

def main():
    parser = argparse.ArgumentParser(description="Run B-cell epitope prediction tools.")
    
    # Add the tool selection argument to accept multiple tools
    parser.add_argument(
        "--tool", 
        choices=["bepipred3", "discotope3", "epigraph"],
        required=True,
        nargs='+',  # Allow multiple tools to be selected
        help="Select one or more tools to run (bepipred3, discotope3, epigraph)"
    )

    # Common argument for temporary output directory
    parser.add_argument(
        "--temp_dir", 
        type=str, 
        default="temp",  # Default temp directory
        help="Temporary directory to store intermediate outputs (default: 'temp')"
    )

    # Common argument for final output directory
    parser.add_argument(
        "--out_dir", 
        type=str, 
        default="output",  # Default output directory
        help="Final output directory (default: 'output')"
    )
    
    # Arguments for Discotope3 & EpiGraph
    parser.add_argument(
        "--pdb", 
        type=str, 
        help="PDB file for Discotope3 & EpiGraph (required)"
    )
    
    # Arguments for Discotope3
    parser.add_argument(
        "--pdb_path", 
        type=str, 
        help="directory where you where you parse the custom pdb and save the pdb downloaded from rcsb.org"
    )

    # Arguments for Bepipred3
    parser.add_argument(
        "--fasta", 
        type=str, 
        help="FASTA file for Bepipred3 (required)"
    )

    parser.add_argument(
        "--pred", 
        type=str, 
        default="vt_pred", 
        help="Prediction model for Bepipred3 (default: 'vt_pred')"
    )

    args = parser.parse_args()

    # Loop through selected tools and run the corresponding function
    for tool in args.tool:
        if tool == "discotope3":
            if not args.pdb:
                print("Error: --pdb is required for Discotope3.")
                sys.exit(1)
            elif not args.pdb_path:
                print("Error: --pdb_path is required for Discotope3.")
                sys.exit(1)
            print(f"Running {tool}...")
            run_discotope3(args.pdb, args.pdb_path, args.temp_dir)
        elif tool == "epigraph":
            if not args.pdb:
                print("Error: --pdb is required for EpiGraph.")
                sys.exit(1)
            elif not args.pdb_path:
                print("Error: --pdb_path is required for EpiGraph.")
                sys.exit(1)
            print(f"Running {tool}...")
            run_epigraph(args.pdb, args.pdb_path, "cuda", args.temp_dir, args.temp_dir)
        elif tool == "bepipred3":
            if not args.fasta:
                print("Error: --fasta is required for Bepipred3.")
                sys.exit(1)
            print(f"Running {tool}...")
            run_bepipred3(args.fasta, args.temp_dir, args.pred)

    # After running tools, standardize the outputs
    print("Standardizing outputs...")
    standardize_outputs(args.temp_dir, args.out_dir)

if __name__ == "__main__":
    main()
