import pytest
import numpy as np
from molbox.rdkit_3d_ops import (
    create_3d_mol_from_smiles,
    get_center_of_mass,
    rotate_conformer,
    set_conformer_positions_,
    set_mol_positions_,
    translate_conformer,
    translate_mol,
)
from scipy.spatial.transform import Rotation


@pytest.fixture
def o2():
    o2 = create_3d_mol_from_smiles("O=O")
    set_mol_positions_(o2, [[0, 0, 0], [1, 1, 1]])
    return o2


@pytest.fixture
def h2o():
    h2o = create_3d_mol_from_smiles("O")
    set_mol_positions_(h2o, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    return h2o


def test_get_center_of_mass(o2):
    center_of_mass = get_center_of_mass(o2)
    assert np.allclose(center_of_mass, [0.5, 0.5, 0.5])


def test_set_conformer_positions_(o2):
    conformer = o2.GetConformer()
    positions = [[1, 1, 1], [2, 2, 2]]
    set_conformer_positions_(conformer, positions)
    assert np.allclose(conformer.GetPositions(), positions)


def test_translate_conformer(o2):
    conformer = translate_conformer(o2.GetConformer(), [1, 1, 1])
    assert np.allclose(conformer.GetPositions(), [[1, 1, 1], [2, 2, 2]])


def test_translate_mol(o2):
    mol = translate_mol(o2, [1, 1, 1])
    assert np.allclose(mol.GetConformer().GetPositions(), [[1, 1, 1], [2, 2, 2]])


def test_rotate_conformer(h2o):
    rotation = Rotation.from_euler("x", 90, degrees=True)
    center = np.array([0, 0, 0])
    conformer = rotate_conformer(h2o.GetConformer(), rotation, center)
    assert np.allclose(conformer.GetPositions(), [[1, 0, 0], [0, 0, 1], [0, -1, 0]])
