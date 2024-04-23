from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation


def create_3d_mol_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    AllChem.EmbedMolecule(mol, params)
    return mol


def get_center_of_mass(mol: Chem.Mol) -> np.ndarray:
    positions = mol.GetConformer().GetPositions()
    masses = np.array([atom.GetMass() for atom in mol.GetAtoms()], dtype=float)
    weighted_positions = masses[:, np.newaxis] * positions
    return weighted_positions.sum(axis=0) / masses.sum()


def set_conformer_positions_(conformer: Chem.Conformer, positions: Any) -> None:
    positions = np.asarray(positions, dtype=float)
    for i in range(conformer.GetNumAtoms()):
        conformer.SetAtomPosition(i, positions[i])


def set_mol_positions_(mol: Chem.Mol, positions: Any) -> None:
    positions = np.asarray(positions, dtype=float)
    set_conformer_positions_(mol.GetConformer(), positions)


def translate_conformer(conformer: Chem.Conformer, translation: Any) -> Chem.Conformer:
    new_conformer = Chem.Conformer(conformer)
    translation = np.asarray(translation, dtype=float)
    positions = new_conformer.GetPositions() + translation
    set_conformer_positions_(new_conformer, positions)
    return new_conformer


def translate_mol(mol: Chem.Mol, translation: Any) -> Chem.Mol:
    new_mol = Chem.Mol(mol)
    conformer = translate_conformer(new_mol.GetConformer(), translation)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(conformer)
    return new_mol


def rotate_conformer(
    conformer: Chem.Conformer, rotation: Rotation, center: np.ndarray
) -> Chem.Conformer:
    new_conformer = Chem.Conformer(conformer)
    new_conformer = translate_conformer(new_conformer, -center)
    positions = rotation.apply(new_conformer.GetPositions())
    set_conformer_positions_(new_conformer, positions)
    new_conformer = translate_conformer(new_conformer, center)
    return new_conformer


def rotate_mol(mol: Chem.Mol, rotation: Rotation, center: Any = None) -> Chem.Mol:
    new_mol = Chem.Mol(mol)
    if center is None:
        center = np.array([0, 0, 0], dtype=float)
    elif center == "COM":
        center = get_center_of_mass(new_mol)
    else:
        center = np.asarray(center, dtype=float)
    conformer = rotate_conformer(new_mol.GetConformer(), rotation, center)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(conformer)
    return new_mol
