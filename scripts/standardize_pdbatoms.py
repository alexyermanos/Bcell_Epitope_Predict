import os
import argparse
import pandas as pd
from Bio import PDB, SeqIO
from Bio import pairwise2

def standardize_pdb(pdb_file, output_dir):
    # Ensure output_dir exists
    os.makedirs(output_dir, exist_ok=True)

    # Compute output CSV path: {output_dir}/{pdb_basename}.csv
    pdb_basename = os.path.splitext(os.path.basename(pdb_file))[0]
    output_csv = os.path.join(output_dir, f"{pdb_basename}_ca_atoms.csv")

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
                        'structure' : pdb_basename,
                        'model': model,
                        'chain': chain.id,
                        'resid': res.get_resname(),
                        'resno': res.id[1],
                        'icode': res.id[2],
                        'atom_name': ca_atom.get_name(),
                        'altloc': ca_atom.get_altloc(),
                        'x': ca_atom.coord[0],
                        'y': ca_atom.coord[1],
                        'z': ca_atom.coord[2],
                        'occupancy': ca_atom.get_occupancy(),
                        'bfactor': ca_atom.get_bfactor(),
                        'element': ca_atom.element
                    })
    atom_df = pd.DataFrame(atom_records)
    print(atom_df)

    # Extract SEQRES sequences for each chain
    seqres_records = list(SeqIO.parse(pdb_file, "pdb-seqres"))
    seqres_data = []
    for rec in seqres_records:
        chain_id = rec.id.split(":")[-1]
        for i, aa in enumerate(str(rec.seq), 1):
            seqres_data.append({
                'structure' : pdb_basename,
                'chain': chain_id,
                'seqresid': aa,
                'seqresno': i
            })
    seqres_df = pd.DataFrame(seqres_data)
    print(seqres_df)

    # Align CA atoms to SEQRES sequences for each chain
    standardized_records = []
    for chain in seqres_df['chain'].unique():
        chain_seqres = seqres_df[seqres_df['chain'] == chain]
        chain_atoms = atom_df[atom_df['chain'] == chain]

        seqres_seq = "".join(chain_seqres['seqresid'])
        print(f'seq: {seqres_seq}')
        
        atom_seq = "".join([PDB.Polypeptide.three_to_one(r) for r in chain_atoms['resid']])
        print(f'atom seq: {atom_seq}')
        
        match_score = 1
        mismatch_score = -1
        gap_open = -1       # minor penalty for opening a gap
        gap_extend = 0      # no extra penalty for extending a gap

        alignments = pairwise2.align.globalms(
            seqres_seq,
            atom_seq,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            penalize_end_gaps=(False, False)  # free gaps at start/end of both sequences
        )
        print(f'alignment: {alignments}')
        
        aln_seqres, aln_atoms, _, _, _ = alignments[0]
        print(f'aln_seqres: {aln_seqres}')
        print(f'aln_atoms: {aln_atoms}')

        seqres_index = 0
        atom_index = 0
        for s, a in zip(aln_seqres, aln_atoms):
            if a != '-':
                atom_row = chain_atoms.iloc[atom_index].to_dict()
                atom_row['seqresid'] = s
                atom_row['seqresno'] = seqres_index + 1
                standardized_records.append(atom_row)
                atom_index += 1
            else:
                standardized_records.append({
                    'structure' : pdb_basename,
                    'chain': chain,
                    'resid': '-',
                    'resno': '-',
                    'seqresid': s,
                    'seqresno': seqres_index + 1
                })
            seqres_index += 1

    # Save CSV
    out_df = pd.DataFrame(standardized_records)
    print(out_df)
    
    out_df.to_csv(output_csv, index=False)
    print(f"Saved standardized atom records to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standardize PDB SEQRES vs CA atoms")
    parser.add_argument("--pdb_file", required=True, help="Input PDB file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    args = parser.parse_args()

    standardize_pdb(args.pdb_file, args.output_dir)
