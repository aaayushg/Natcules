from pathlib import Path

import pandas as pd
from deepchem.data import NumpyDataset
from deepchem.feat import RDKitDescriptors
from joblib import load
from rdkit import Chem

# === Config ===
INPUT_CSV = "data/enumerated_smiles.csv"
OUTPUT_CSV = "results/solubility_predictions.csv"
SMILES_COLUMN = "smiles"

# === Load SMILES ===
df = pd.read_csv(INPUT_CSV)
df = df[df[SMILES_COLUMN].notna()]
valid_smiles = [s for s in df[SMILES_COLUMN] if Chem.MolFromSmiles(s)]
df = df[df[SMILES_COLUMN].isin(valid_smiles)].reset_index(drop=True)

# === Featurize ===
featurizer = RDKitDescriptors()
features = featurizer.featurize(df[SMILES_COLUMN].tolist())
dataset = NumpyDataset(X=features)

# === Load pretrained model (or substitute with your own) ===
# Example: load a local .joblib model

model = load("models/solubility/solubility_rf_model.joblib")
preds = model.predict(features)

# === Save ===
df["predicted_solubility"] = preds
Path("results").mkdir(exist_ok=True)
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved solubility predictions to {OUTPUT_CSV}")

