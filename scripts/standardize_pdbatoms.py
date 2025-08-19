import argparse
import pandas as pd
from Bio import PDB, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

def standardize_pdb(pdb_file, output_csv):
    # Load PDB structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    # Extract CA atoms for each chain
    atom_records = []
    for model in structure:
        for chain in model:
            for res in chain:
                if PDB.is_aa(res, standard=True) and 'CA' in res:
                    ca_atom = res['CA']
                    atom_records.append({
                        'chain': chain.id,
                        'resid': res.get_resname(),
                        'resno': res.id[1]
                    })

    atom_df = pd.DataFrame(atom_records)
    
    # Extract SEQRES sequences for each chain
    seqres_records = list(SeqIO.parse(pdb_file, "pdb-seqres"))
    seqres_data = []
    for rec in seqres_records:
        chain_id = rec.id
        for i, aa in enumerate(str(rec.seq), 1):
            seqres_data.append({
                'chain': chain_id,
                'seqresid': aa,
                'seqresno': i
            })
    seqres_df = pd.DataFrame(seqres_data)

    # Align CA atoms to SEQRES sequences for each chain
    standardized_records = []
    for chain in seqres_df['chain'].unique():
        chain_seqres = seqres_df[seqres_df['chain'] == chain]
        chain_atoms = atom_df[atom_df['chain'] == chain]

        # Create sequences
        seqres_seq = "".join(chain_seqres['seqresid'])
        atom_seq = "".join([PDB.Polypeptide.three_to_one(r) for r in chain_atoms['resid']])

        # Align with no end-gap penalty
        alignments = pairwise2.align.globalms(seqres_seq, atom_seq, 1, -1, 0, 0)
        best_alignment = alignments[0]
        aln_seqres, aln_atoms, _, _, _ = best_alignment

        # Map seqresno/resid
        seqres_index = 0
        atom_index = 0
        for s, a in zip(aln_seqres, aln_atoms):
            if a != '-':
                atom_row = chain_atoms.iloc[atom_index].to_dict()
                atom_row['seqresid'] = s
                atom_row['seqresno'] = seqres_index + 1
                standardized_records.append(atom_row)
                atom_index += 1
            elif a == '-':
                # Missing atom, add placeholder row
                atom_row = {
                    'chain': chain,
                    'resid': '-',
                    'resno': '-',
                    'seqresid': s,
                    'seqresno': seqres_index + 1
                }
                standardized_records.append(atom_row)
            seqres_index += 1

    # Save CSV
    out_df = pd.DataFrame(standardized_records)
    out_df.to_csv(output_csv, index=False)
    print(f"Saved standardized atom records to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standardize PDB SEQRES vs CA atoms")
    parser.add_argument("--pdb_file", required=True, help="Input PDB file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    args = parser.parse_args()

    standardize_pdb(args.pdb_file, args.output_csv)
