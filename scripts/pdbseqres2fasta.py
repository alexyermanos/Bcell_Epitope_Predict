import os
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Extract SEQRES from PDB to FASTA")
    parser.add_argument("pdb_file", type=str, help="Input PDB file")
    parser.add_argument("fasta_file", type=str, help="Output FASTA file")
    args = parser.parse_args()

    # Try SEQRES first
    records = list(SeqIO.parse(args.pdb_file, "pdb-seqres"))

    # Fall back to ATOM if no SEQRES found
    if not records:
        print(f"Warning: No SEQRES records found in {args.pdb_file}, falling back to ATOM records.", flush=True)
        records = list(SeqIO.parse(args.pdb_file, "pdb-atom"))

    if not records:
        raise ValueError(f"No sequence could be extracted from {args.pdb_file}")

    # Write FASTA, enforcing header = PDBID:CHAIN
    with open(args.fasta_file, "w") as out_fasta:
        for rec in records:
            header = rec.id.split()[0]   # keep only "3CX3:A"
            out_fasta.write(f">{header}\n{str(rec.seq)}\n")

    print(f"Saved FASTA to {args.fasta_file}", flush=True)

if __name__ == "__main__":
    main()
