from typing import Optional

from rdkit import Chem

from molbox.box import Box


class MolBox:
    def __init__(
        self, mol: Optional[Chem.Mol] = None, box: Optional[Box] = None
    ) -> None:
        self.mol = Chem.Mol() if mol is None else mol
        self.box = self.get_minimal_box() if box is None else box

    def get_minimal_box(self) -> Box:
        positions = self.mol.GetConformer().GetPositions()
        corner_lo = positions.min(axis=0)
        corner_hi = positions.max(axis=0)
        return Box(corner_lo, corner_hi)

    def add_mol(self, mol: Chem.Mol) -> None:
        self.mol = Chem.CombineMols(self.mol, mol)

    def merge_molbox(self, other: "MolBox") -> None:
        self.add_mol(other.mol)
