import os
import time

import pandas as pd
import requests
from rdkit import Chem

# Folder to save .mol2 files
output_dir = "mol2_files"
os.makedirs(output_dir, exist_ok=True)

base_url = "https://www.tcmsp-e.com/tcmspmol/MOL{:06d}.mol2"
headers = {
    "User-Agent": "Mozilla/5.0"
}

for i in range(15001):  # From MOL000000 to MOL015000
    mol_id = f"{i:06d}"
    url = base_url.format(i)
    output_path = os.path.join(output_dir, f"MOL{mol_id}.mol2")

    try:
        response = requests.get(url, headers=headers, timeout=10)
        if response.status_code == 200 and response.text.strip():
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(response.text)
            print(f"[{mol_id}] Downloaded")
        else:
            print(f"[{mol_id}] Not found or empty")
    except Exception as e:
        print(f"[{mol_id}] Error: {e}")
    time.sleep(0.2)  # Throttle to be respectful to the server


# Directory containing .mol2 files
mol2_dir = "mol2_files"
output = []

for filename in os.listdir(mol2_dir):
    if filename.endswith(".mol2"):
        file_path = os.path.join(mol2_dir, filename)
        try:
            mol = Chem.MolFromMol2File(file_path, sanitize=True)
            if mol:
                smiles = Chem.MolToSmiles(mol)
                output.append({"filename": filename, "smiles": smiles})
            else:
                print(f"Could not parse {filename}")
        except Exception as e:
            print(f"Error processing {filename}: {e}")

# Save to CSV
df = pd.DataFrame(output)
df.to_csv("tcm_smiles.csv", index=False)


