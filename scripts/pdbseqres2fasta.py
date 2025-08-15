import os
from Bio import SeqIO

import argparse

parser = argparse.ArgumentParser(description="Extract SEQRES from PDB to FASTA")
parser.add_argument("pdb_file", type=str, help="Input PDB file")
parser.add_argument("fasta_file", type=str, help="Output FASTA file")
args = parser.parse_args()

# Parse PDB SEQRES
records = list(SeqIO.parse(args.pdb_file, "pdb-seqres"))

# Write FASTA with only record.id as header
with open(args.fasta_file, "w") as out_fasta:
    for rec in records:
        out_fasta.write(f">{rec.id}\n{str(rec.seq)}\n")

print(f"Saved FASTA to {args.fasta_file}")