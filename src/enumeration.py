from rdkit import Chem
from rdkit.Chem import BRICS

def smiles_to_mol(smiles):
    """
    Convert a SMILES string to an RDKit Mol object.
    """
    mol = Chem.MolFromSmiles(smiles)
    return mol

def mol_to_smiles(mol):
    """
    Convert RDKit Mol object to canonical SMILES string.
    """
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)

def enumerate_molecules_from_fragments(fragments, max_output=20):
    """
    Generate new molecules by combining BRICS fragments.
    
    Parameters:
        fragments (list of str): List of fragment SMILES
        max_output (int): Maximum number of molecules to return

    Returns:
        List of SMILES strings representing new molecules
    """
    # Convert SMILES to mols
    mols = [Chem.MolFromSmiles(f) for f in fragments if Chem.MolFromSmiles(f) is not None]
    if not mols:
        return []

    # Build new molecules
    try:
        new_mols = BRICS.BRICSBuild(mols)
        new_smiles = []
        for mol in new_mols:
            if len(new_smiles) >= max_output:
                break
            smi = mol_to_smiles(mol)
            if smi and smi not in new_smiles:
                new_smiles.append(smi)
        return new_smiles
    except Exception as e:
        print(f"Error during enumeration: {e}")
        return []


