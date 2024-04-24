import re

from openbabel import openbabel
from rdkit import Chem

from molbox.molbox import MolBox


def write_lammps_data(molbox: MolBox, filename: str) -> None:
    mol_block = Chem.MolToMolBlock(molbox.mol)

    ob_mol = openbabel.OBMol()
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("mol", "lmpdat")
    ob_conversion.ReadString(ob_mol, mol_block)
    lmpdata_block = ob_conversion.WriteString(ob_mol)

    corner_lo = molbox.box.corner_lo
    corner_hi = molbox.box.corner_hi
    box_section = (
        f"  {corner_lo[0]: .5f}   {corner_hi[0]: .5f} xlo xhi\n"
        f"  {corner_lo[1]: .5f}   {corner_hi[1]: .5f} ylo yhi\n"
        f"  {corner_lo[2]: .5f}   {corner_hi[2]: .5f} zlo zhi"
    )
    pattern = r".*xlo xhi\n.*ylo yhi\n.*zlo zhi"
    lmpdata_block = re.sub(pattern, box_section, lmpdata_block)

    with open(filename, "w") as f:
        f.write(lmpdata_block)
