from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import Recap
from rdkit.Chem import Draw

def smiles_to_mol(smiles):
    """
    Convert SMILES to RDKit Mol object.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    return mol

def brics_fragmentation(smiles):
    """
    Fragment molecule using BRICS method.
    Returns a list of fragment SMILES.
    """
    mol = smiles_to_mol(smiles)
    if mol is None:
        return []
    fragments = BRICS.BRICSDecompose(mol)
    return list(fragments)

def recap_fragmentation(smiles):
    """
    Fragment molecule using RECAP method.
    Returns a list of fragment SMILES.
    """
    mol = smiles_to_mol(smiles)
    if mol is None:
        return []
    recap_tree = Recap.RecapDecompose(mol)
    return list(recap_tree.GetLeaves().keys())

def fragment_molecule(smiles, method='brics'):
    """
    General interface to fragment a molecule.
    method: 'brics' or 'recap'
    """
    if method == 'brics':
        return brics_fragmentation(smiles)
    elif method == 'recap':
        return recap_fragmentation(smiles)
    else:
        raise ValueError("Unsupported fragmentation method. Choose 'brics' or 'recap'.")

# from fragmentation import fragment_molecule

# smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
# fragments = fragment_molecule(smiles, method='brics')
# print(fragments)