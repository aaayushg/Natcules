import pandas as pd
import torch
from rdkit import Chem
from torch.utils.data import DataLoader
from pathlib import Path

# Import DeepTox modules (assumes cloned as 'DeepTox/')
import sys
sys.path.append("DeepTox")
from utils import load_model, get_test_transforms, collate_fn
from dataset import Tox21

# === CONFIG ===
MODEL_PATH = "DeepTox/models/best_model.pth"
INPUT_CSV = "data/enumerated_smiles.csv"  # Your input file
OUTPUT_CSV = "results/toxicity_predictions.csv"
SMILES_COLUMN = "smiles"
BATCH_SIZE = 32

# === LOAD DATA ===
df = pd.read_csv(INPUT_CSV)
smiles_list = df[SMILES_COLUMN].dropna().tolist()

# RDKit filter
valid_smiles = [s for s in smiles_list if Chem.MolFromSmiles(s) is not None]
df = df[df[SMILES_COLUMN].isin(valid_smiles)].reset_index(drop=True)

# === Prepare DeepTox model ===
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = load_model(MODEL_PATH, device=device)
model.eval()

# === DeepTox Dataset ===
dataset = Tox21(smiles=df[SMILES_COLUMN].tolist(), labels=None, transform=get_test_transforms())
dataloader = DataLoader(dataset, batch_size=BATCH_SIZE, collate_fn=collate_fn)

# === Run predictions ===
all_preds = []
with torch.no_grad():
    for batch in dataloader:
        inputs = batch["features"].to(device)
        logits = model(inputs)
        preds = torch.sigmoid(logits).cpu().numpy()
        all_preds.extend(preds)

# Add predictions to DataFrame
endpoint_labels = [f"Tox21_{i}" for i in range(len(all_preds[0]))]
pred_df = pd.DataFrame(all_preds, columns=endpoint_labels)
result_df = pd.concat([df, pred_df], axis=1)

# Save results
Path("results").mkdir(exist_ok=True)
result_df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved toxicity predictions to {OUTPUT_CSV}")

