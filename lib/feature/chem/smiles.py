import rdkit import Chem

def molfile_to_smiles(mol_file):
    mol = Chem.MolFromMolFile(mol_file)
    if mol is not None:
        smiles = Chem.MolToSmiles(mol)
        return smiles
    else:
        return "Error: Unable to parse mol file."
