from rdkit import Chem
from rdkit.Chem import Descriptors


def calculate_descriptors(smiles):
    """
    Calculate common RDKit-based molecular descriptors.
    Returns a dictionary of descriptor values.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    descriptors = {
        'MolWt': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'RingCount': Descriptors.RingCount(mol)
    }
    return descriptors

