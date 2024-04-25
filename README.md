# molbox

Molbox is a simple Python package that extends [rdkit](https://www.rdkit.org/) with spatial [transformations](src/molbox/rdkit_3d_ops.py) of 3D molecules and semantics for multiple molecules in a rectangular [box](src/molbox/molbox.py). It allows you to quickly create simulation boxes and write them to data files that can be read by [LAMMPS](https://www.lammps.org/).

## Installation

```
pip install molbox
```

## Usage

```python
from molbox import Box, MolBox
from molbox.lammps_data import write_lammps_data
from molbox.rdkit_3d_ops import create_3d_mol_from_smiles, rotate_mol
from rdkit import Chem
from scipy.spatial.transform import Rotation


def main():
    box = Box([0, 0, 0], [175, 175, 175])
    o2_box = MolBox(box=box)
    o2_mol = create_3d_mol_from_smiles("O=O")
    num_mols = 128
    for _ in range(num_mols):
        rotation = Rotation.random()
        new_mol = rotate_mol(o2_mol, rotation, "COM")
        o2_box.add_mol_randomly(new_mol)
    Chem.MolToMolFile(o2_box.mol, "examples/o2_box.mol")
    write_lammps_data(o2_box, "examples/o2_box.lmpdat")


if __name__ == "__main__":
    main()
```
