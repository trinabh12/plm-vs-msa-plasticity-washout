import os
import requests
import urllib.request
import numpy as np
import time
from Bio import PDB

# Output directory relative to src/
PREDICTED_DIR = "../data/processed/predicted_structures/"
ESMFOLD_DIR = os.path.join(PREDICTED_DIR, "esmfold")
ALPHAFOLD_DIR = os.path.join(PREDICTED_DIR, "alphafold")

os.makedirs(ESMFOLD_DIR, exist_ok=True)
os.makedirs(ALPHAFOLD_DIR, exist_ok=True)

ESMATLAS_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

def predict_esmfold(protein_id: str, sequence: str) -> str:
    """Queries the ESMFold API with a retry mechanism for 504 timeouts."""
    # Dropped to 150aa to bypass Meta's current server timeouts
    query_seq = sequence[:150] 
    out_path = os.path.join(ESMFOLD_DIR, f"{protein_id}_esmfold.pdb")
    
    if os.path.exists(out_path):
        print(f"[CACHE] Found existing ESMFold structure: {out_path}")
        return out_path

    # Retry loop to handle 504 Gateway Timeouts gracefully
    for attempt in range(3):
        print(f"Predicting ESMFold structure (Attempt {attempt+1}) for {protein_id} (Length: {len(query_seq)} aa)...")
        response = requests.post(ESMATLAS_URL, data=query_seq)
        
        if response.status_code == 200:
            with open(out_path, "w") as f:
                f.write(response.text)
            print(f"SUCCESS: Saved ESMFold structure -> {out_path}")
            return out_path
        elif response.status_code == 504:
            print("Server overloaded (504). Waiting 3 seconds before retrying...")
            time.sleep(3)
        else:
            print(f"ERROR: ESMFold API returned status code {response.status_code}")
            break
            
    return None

def fetch_alphafold_db(uniprot_id: str) -> str:
    """Retrieves structure dynamically via the AlphaFold API to prevent 404s."""
    out_path = os.path.join(ALPHAFOLD_DIR, f"AF_{uniprot_id}.pdb")
    if os.path.exists(out_path):
        print(f"[CACHE] Found existing AlphaFold structure: {out_path}")
        return out_path
        
    print(f"Querying AlphaFold API for UniProt ID {uniprot_id}...")
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    
    try:
        response = requests.get(api_url)
        if response.status_code == 200:
            data = response.json()
            pdb_url = data[0]['pdbUrl'] # Grabs the exact URL for Fragment 1
            print(f"Resolved dynamic URL: {pdb_url}")
            
            urllib.request.urlretrieve(pdb_url, out_path)
            print(f"SUCCESS: Saved AlphaFold structure -> {out_path}")
            return out_path
        else:
            print(f"ERROR: AlphaFold API returned {response.status_code}")
            return None
    except Exception as e:
        print(f"ERROR fetching AlphaFold DB entry for {uniprot_id}: {e}")
        return None

def extract_plddt(pdb_path: str) -> list[float]:
    """Parses B-factor column (pLDDT confidence) for C-alpha atoms in a PDB."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    return [atom.get_bfactor() for atom in structure.get_atoms() if atom.get_name() == "CA"]

def calculate_ca_rmsd(pdb_path1: str, pdb_path2: str) -> float:
    """Aligns two structures in 3D space and calculates the true RMSD."""
    parser = PDB.PDBParser(QUIET=True)
    s1 = parser.get_structure("struct1", pdb_path1)
    s2 = parser.get_structure("struct2", pdb_path2)
    
    # Extract Alpha Carbon atoms
    ca1 = [a for a in s1.get_atoms() if a.get_name() == "CA"]
    ca2 = [a for a in s2.get_atoms() if a.get_name() == "CA"]
    
    # Match lengths for comparison (using the 150aa cutoff)
    min_len = min(len(ca1), len(ca2))
    ca1 = ca1[:min_len]
    ca2 = ca2[:min_len]
    
    # Superimpose (align) Structure 2 onto Structure 1
    super_imposer = PDB.Superimposer()
    super_imposer.set_atoms(ca1, ca2)
    super_imposer.apply(s2.get_atoms())
    
    # Return the aligned RMSD
    return float(super_imposer.rms)