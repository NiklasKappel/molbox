import numpy as np
import pytest
from molbox.rdkit_helpers import create_3d_mol_from_smiles, set_conformer_positions


@pytest.fixture
def mol():
    return create_3d_mol_from_smiles("O=O")


def test_set_conformer_positions(mol):
    conformer = mol.GetConformer()
    positions = [[0, 0, 0], [1, 1, 1]]
    set_conformer_positions(conformer, positions)
    assert np.allclose(conformer.GetPositions(), positions)
