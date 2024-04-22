from molbox import Box, MolBox
from molbox.rdkit_helpers import create_3d_mol_from_smiles
from rdkit import Chem


def main():
    box = Box([0, 0, 0], [10, 10, 10])
    o2_box = MolBox(box=box)
    o2_mol = create_3d_mol_from_smiles("O=O")
    num_mols = 32
    for _ in range(num_mols):
        new_mol = Chem.Mol(o2_mol)
        # rotate the molecule randomly
        # translate the molecule randomly within the box
        o2_box.add_mol(new_mol)


if __name__ == "__main__":
    main()
