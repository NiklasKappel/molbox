import copy
from typing import Optional

import numpy as np
from rdkit import Chem

from molbox.box import Box
from molbox.rdkit_3d_ops import translate_mol


class MolBox:
    def __init__(
        self, mol: Optional[Chem.Mol] = None, box: Optional[Box] = None
    ) -> None:
        self.mol = Chem.Mol() if mol is None else Chem.Mol(mol)
        self.box = self.get_minimal_box() if box is None else copy.copy(box)

    def get_minimal_box(self) -> Box:
        positions = self.mol.GetConformer().GetPositions()
        corner_lo = positions.min(axis=0)
        corner_hi = positions.max(axis=0)
        return Box(corner_lo, corner_hi)

    def add_mol(self, mol: Chem.Mol) -> None:
        self.mol = Chem.CombineMols(self.mol, mol)

    def add_mol_randomly(self, mol: Chem.Mol) -> None:
        rand = np.random.rand(3)
        other = MolBox(mol)
        # fmt: off
        translation = (
            (1 - rand) * (self.box.corner_lo - other.box.corner_lo)
            + rand * (self.box.corner_hi - other.box.corner_hi)
        )
        # fmt: on
        self.add_mol(translate_mol(mol, translation))
