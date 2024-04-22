from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def create_3d_mol_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    AllChem.EmbedMolecule(mol, params)
    return mol


def set_conformer_positions(conformer: Chem.Conformer, positions: Any) -> None:
    positions = np.asarray(positions, dtype=float)
    for i in range(conformer.GetNumAtoms()):
        conformer.SetAtomPosition(i, positions[i])


def set_mol_positions(mol: Chem.Mol, positions: np.ndarray) -> None:
    conformer = mol.GetConformer()
    set_conformer_positions(conformer, positions)
