import argparse
import pandas as pd
import numpy as np
import shutil
import wget
import os
import sys
import warnings

epigraph_path = os.path.join(os.path.dirname(__file__), "EpiGraph")
sys.path.insert(0, epigraph_path)
sys.dont_write_bytecode = True

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.DSSP import dssp_dict_from_pdb_file, residue_max_acc
from torch_geometric.data import Data
from esm_embedding import esm_if_2_embedding

from epigraph_model import GAT
import torch

warnings.filterwarnings('ignore')

arg_parser = argparse.ArgumentParser(description="Data fetch and Processing")

arg_parser.add_argument("--pdb", type=str, help="query PDB", required=True)

arg_parser.add_argument("--pdb_path", type=str, default="PDB", help="directory where you where you parse the custom pdb and save the pdb downloaded from rcsb.org")

arg_parser.add_argument("--save_path", type=str, default="PDB_Processed", help="directory where you save the heteroatom removed pdb")

arg_parser.add_argument("--distance_threshold", type=int, default=10, help="distance_threshold that generates edge connection between nodes")

arg_parser.add_argument("--RSA_threshold", type=float, default=0.15, help="nodes above RSA_threshold are train_mask==True, below are train_mask==False")

arg_parser.add_argument("--model_path", type=str, default="checkpoint", help="directory where models are saved ")

arg_parser.add_argument("--kfold", type=int, default=10, help="number of models to ensemble")

arg_parser.add_argument("--out_path", type=str, default="Result", help="directory where the pred.csv is saved")

arg_parser.add_argument("--classification_threshold", type=float, default=0.1481, help="threshold for classification")

arg_parser.add_argument("--device", type=str, default="cuda")

args = arg_parser.parse_args()

pdb = args.pdb
#pdb = pdb.lower()
pdb_path = args.pdb_path
save_path = args.save_path
distance_threshold = args.distance_threshold
RSA_threshold = args.RSA_threshold
csv_out = args.pdb
model_path = args.model_path
kfold = args.kfold
out_path = args.out_path
classification_threshold = args.classification_threshold
device = torch.device(args.device if torch.cuda.is_available() else "cpu")


#if os.path.isdir(pdb_path):
#    pass
#else:
#    os.mkdir(pdb_path)

os.makedirs(pdb_path, exist_ok=True) #NEW

pdb_file_path = os.path.join(pdb_path, f"{pdb}.pdb")

if os.path.exists(pdb_file_path):
    print(f"Found local PDB file: {pdb_file_path}")
else:
    print(f"{pdb_file_path} not found locally. Downloading from rcsb.org...")
    print()
    try:
        pdb_lower = pdb.lower()
        wget.download(f"https://files.rcsb.org/download/{pdb_lower}.pdb")
        shutil.copy(f"{pdb_lower}.pdb", os.path.join(pdb_path, f"{pdb}.pdb"))
        os.remove(f"{pdb_lower}.pdb")
        print(f"Downloaded and saved to {pdb_file_path}")
    except Exception as e:
        print("="*50)
        print("Error occurred.", e)
        print(f"{pdb_lower}.pdb not found in rcsb.org and not found locally")
        print("Please check the query is available in rcsb.org as pdb format")
        print("="*50)
        sys.exit(1)
    
# class that removes hetero atoms in PDB format
# ref https://stackoverflow.com/questions/25718201/remove-heteroatoms-from-pdb

print("\nData Processing...")
class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

# save pre-processed pdb format into save_path
try:
    Bio_parser = PDBParser()
    model = Bio_parser.get_structure(f"{pdb}", f"{pdb_path}/{pdb}.pdb")
except Exception as e:
    print("="*50)
    print("Error occured.", e)
    print(f"{pdb}.pdb not found in {pdb_path}/")
    print(f"Please check the query PDB is available in {pdb_path} as pdb format")
    print("="*50)
    pass

if os.path.isdir(save_path):
    pass
else:
    os.mkdir(save_path)

io = PDBIO()
io.set_structure(model)
io.save(f"{save_path}/{pdb}.pdb", NonHetSelect())

def euclidean_dist(x, y):
    return ((x[:, None] - y) ** 2).sum(-1).sqrt()

def edge_connection(coord_list, threshold):
    # Compute pairwise euclidean distances
    distances = euclidean_dist(coord_list, coord_list)
    
    # to avoid self-connection, make the distance 0 between self nodes into infinity
    distances.fill_diagonal_(float("inf"))

    # edges are constructed within threshold 
    edges = (distances < threshold).nonzero(as_tuple=False).t()
    
    return edges

def generate_graph(pdb, save_path, distance_threshold, RSA_threshold):    
    esm_if_rep, esm2_rep, node_list, coord_list = esm_if_2_embedding(pdb, save_path)
    
    esm_node_features = torch.concat((esm_if_rep, esm2_rep), dim=1)

    node_all_list = []
    for chain_node in node_list:
        for node in chain_node:
            node_all_list.append(node)

    coord_all_list = []
    for chain_coord in coord_list:
        for coord in chain_coord:
            coord_all_list.append(coord)

    coord_all_list = torch.tensor(np.array(coord_all_list))

    edges = edge_connection(coord_all_list, threshold=distance_threshold)

    dssp = dssp_dict_from_pdb_file(f"{save_path}/{pdb}.pdb", DSSP="mkdssp")
    
#    print("Available DSSP keys:")
#    for key in dssp[0].keys():
#        print(key)


    rsa_list = []
    for node in node_all_list:
        chain, res_name, res_id = node.split(":")
        try:
            # indexing the dssp such as ('A', (' ', 53, ' '))
            key = (chain, (' ', int(res_id), ' '))
            
            # generate rsa by normalizing asa by residue_max_acc -> 
            rsa = dssp[0][key][2] / residue_max_acc["Sander"][res_name]
            rsa_list.append(rsa)
        except:
            print("Key Error... appending rsa: 0")
            rsa_list.append(0)
        
        # The surface residues were selected with certain RSA cutoff 
        # surface residues above RSA cutoff is True, buried residues below RSA cutoff is False


    train_mask = torch.tensor([rsa >= RSA_threshold for rsa in rsa_list])

    data = Data(
        coords=coord_all_list,
        node_id=node_all_list,
        x=esm_node_features,
        edge_index=edges.contiguous(),
        num_nodes=len(node_all_list),
        name=pdb,
        train_mask=train_mask,
        rsa=rsa_list
    )

    if data.x.shape[0] != data.num_nodes:
        print("=" * 50)
        print(f"pdb {data.name} got an error; node assignment error")
        print("=" * 50)

    return data

data = generate_graph(pdb, save_path, distance_threshold=distance_threshold, RSA_threshold=RSA_threshold)

# model inference with ensemble model
os.makedirs(model_path, exist_ok=True)

def ensemble_pred(model_path, data, kfold, RSA_threshold, device):
    model = GAT(in_dim=1792, hid_dim=128, out_dim=1, num_head=8, out_head=1)

    model.to(device)
    
    num_nodes = data.num_nodes
    pred_ensem = torch.zeros([num_nodes])
    
    pt_list = [pt for pt in os.listdir(model_path) if pt.endswith(".pt")]
    if not pt_list:
        raise RuntimeError(f"No model (.pt) files found in {model_path}")
        
    #print(pt_list)
        
    for pt in pt_list:
        model.load_state_dict(torch.load(f'{model_path}/{pt}', map_location=device), strict = False)
        model.eval()
        
        rsa_list = []
        node_list = []
        pred_label_list = []
        with torch.no_grad():
            
            data.to(device)
            out = model(data)
            y_pred = out.reshape(-1)
                
            for pred, rsa, node_id in zip(y_pred, data.rsa, data.node_id):
                # residue lower than RSA_threshold; regarded as buried residue; no epitope
                if rsa < RSA_threshold:
                    pred_label_list.append(0)
                else:
                    pred_label_list.append(pred)
                rsa_list.append(rsa)
                node_list.append(node_id)
            
        pred_label_list = torch.tensor(pred_label_list)
        pred_ensem += pred_label_list
        
    pred_ensem_list = []
    for i in (pred_ensem / kfold): # 10 fold
        pred_ensem_list.append(i.cpu().item())

    return pred_ensem_list, rsa_list, node_list

ensem_pred, rsa_list, node_list = ensemble_pred(model_path, data, kfold=kfold, RSA_threshold=RSA_threshold, device=device)

def to_labels(ensem_pred, threshold):
    ensem_binary = []
    for pred in ensem_pred:
        if pred > threshold:
            ensem_binary.append(1)
        else:
            ensem_binary.append(0)
    
    return ensem_binary

ensem_binary = to_labels(ensem_pred, classification_threshold)

# save the inference result in csv file 

data = {"PDB": pdb,
        "Residue": node_list,
        "Score": ensem_pred,
        "Epitope": ensem_binary,
        "RSA": rsa_list}

df = pd.DataFrame(data)

#if os.path.isdir(out_path):
#    pass
#else:
#    os.mkdir(out_path)

os.makedirs(out_path, exist_ok=True) #NEW

df.to_csv(f"{out_path}/{csv_out}.csv")

print(f"{out_path}/{csv_out}.csv saved!")
