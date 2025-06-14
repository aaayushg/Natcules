from rdkit import Chem
from rdkit.Chem import QED
from descriptor_generator import calculate_descriptors

def predict_qed(smiles):
    """
    Predict QED (quantitative estimate of drug-likeness) score.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return QED.qed(mol)

def predict_properties(smiles):
    """
    Return a full dictionary of predicted properties.
    Includes RDKit descriptors and QED.
    """
    desc = calculate_descriptors(smiles)
    if desc is None:
        return None
    desc['QED'] = predict_qed(smiles)

    # Placeholder for additional model predictions (ADME/Tox, etc.)
    # desc['ToxicityScore'] = custom_model.predict(smiles)

    return desc

# from property_predictor import predict_properties

# smiles = "CC(=O)Oc1ccccc1C(=O)O"
# props = predict_properties(smiles)
# print(props)
