#!/usr/bin/env python3

import os
import argparse
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP, dssp_dict_from_pdb_file, residue_max_acc

print(residue_max_acc["Sander"])

def run_dssp(pdb, pdb_path="PDB", dssp_exec="mkdssp", dssp_path="temp"):
    """
    Run DSSP (mkdssp) on a given PDB file and save the output.

    Parameters:
    pdb (str): PDB filename without extension (e.g., '2BIB')
    pdb_path (str): Directory where the PDB file is located
    dssp_exec (str): DSSP executable (default 'mkdssp')
    dssp_path (str): Directory to save DSSP output (default 'temp')
    """
    # Ensure DSSP output directory exists
    os.makedirs(dssp_path, exist_ok=True)

    pdb_file = os.path.join(pdb_path, f"{pdb}.pdb")

    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb, pdb_file)
    model = structure[0]  # Use first model

    # Run DSSP
    dssp = DSSP(model, pdb_file, dssp=dssp_exec)

    # Print a few keys and their data
    print("Example DSSP entries from DSSP:")
    for i, key in enumerate(dssp.keys()):
        if i >= 5:  # print only first 5 residues
            break
        print(f"{key}: {dssp[key]}")
    
    dssp_tuple = dssp_dict_from_pdb_file(pdb_file, DSSP=dssp_exec)
    dssp_dict = dssp_tuple[0]  
    
    print("Example DSSP entries from dssp_dict_from_pdb_file:")
    for i, key in enumerate(dssp_dict.keys()):
        if i >= 5:
          break
        print(f"{key}: {dssp_dict[key]}")
    
    return dssp

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DSSP on a PDB ID and save output")
    parser.add_argument("--pdb", type=str, required=True, help="PDB ID (e.g., 2BIB)")
    parser.add_argument("--pdb_path", type=str, default="PDB", help="Directory with PDB files")
    parser.add_argument("--dssp_exec", type=str, default="mkdssp", help="DSSP executable")
    parser.add_argument("--dssp_path", type=str, default="temp", help="Directory to save DSSP output")

    args = parser.parse_args()

    run_dssp(args.pdb, args.pdb_path, args.dssp_exec, args.dssp_path)
