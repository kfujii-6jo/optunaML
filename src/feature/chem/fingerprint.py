from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import MACCSkeys
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
from rdkit.Chem.AllChem import GetHashedTopologicalTorsionFingerprintAsBitVect
import numpy as np
import pandas as pd

def mol_to_ecfp(molecule, radius=3, nBits=1024):
    ecfp = AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=nBits)
    ecfp_list = [int(bit) for bit in ecfp.ToBitString()]
    return ecfp_list

def mol_to_fcfp(molecule, radius=2, nBits=1024):
    fcfp = AllChem.GetMorganFingerprintAsBitVect(molecule, radius, useFeatures=True, nBits=nBits)
    fcfp_list = [int(bit) for bit in fcfp.ToBitString()]
    return fcfp_list

def mol_to_maccs_keys(molecule):
    maccs_keys = MACCSkeys.GenMACCSKeys(molecule)
    maccs_keys_list = [int(bit) for bit in maccs_keys.ToBitString()]
    return maccs_keys_list

def mol_to_avalon_fp(molecule):
    avalon_fp = GetAvalonFP(molecule)
    avalon_fp_list = [int(bit) for bit in avalon_fp.ToBitString()]
    return avalon_fp_list

def mol_to_topological_fp(molecule):
    topological_fp = GetHashedTopologicalTorsionFingerprintAsBitVect(molecule)
    topological_fp_list = [ 1 if topological_fp.GetBit(i) else 0 for i in range(topological_fp.GetNumBits())]
    return topological_fp_list

def smiles_to_ecfp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol_to_ecfp(mol)

def smiles_to_fcfp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol_to_fcfp(mol)

def smiles_to_maccs_keys(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol_to_maccs_keys(mol)

def smiles_to_avalon(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol_to_avalon_fp(mol)

def smiles_to_topological(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol_to_topological_fp(mol)
