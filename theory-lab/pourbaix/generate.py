"""
Generate molecules and material structures and optimize them using MACE calculator.
"""

import numpy as np
from ase.build import molecule
from ase.cell import Cell
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import UnitCellFilter
from mace.calculators import mace_mp

# Define molecules
molecules = [
    molecule("O2"),
    molecule("H2"),
    molecule("H2O"),
]

for m in molecules:
    m.cell = Cell.new([10, 10, 10])
    m.pbc = True
    m.info["type"] = "molecule"
    if sum(np.abs(m.get_initial_magnetic_moments())) > 0:
        m.info["spin_polarized"] = True
        del m.arrays["initial_magmoms"]
    m.info["name"] = m.get_chemical_formula()

# Read structures from Crystallography Open Database cif files
cifs = {
    "Li": 9008473,
    "Na": 9008545,
    "Li2O": 1010064,
    "Na2O": 1010876,
    "Co_heterogenite_2h": 9009449,
    "Co_heterogenite_3r": 9009884,
    "Co_hydroxide": 1548810,
    "CoO_1528838": 1528838,
    "CoO_1533087": 1533087,
    "CoO2": 3000496,
    "Co3O4": 1538531,
    "Co": 1534891,
}

materials = []

for name, cod_id in cifs.items():
    m = read(name + ".cif", format="cif")
    m.info["type"] = "material"
    m.info["name"] = name
    m.info["cod_id"] = cod_id
    if "occupancy" in m.info:
        del m.info["occupancy"]
    if "spacegroup_kinds" in m.arrays:
        del m.arrays["spacegroup_kinds"]
    materials.append(m)

# Calculate the energies with MACE-MP
structures = molecules + materials
calculator = mace_mp(model="medium", default_dtype="float64")

for s in structures:
    s.calc = calculator
    s.info["MACE_energy"] = s.get_potential_energy()
    print(s.info["MACE_energy"])

write("vasp/structures.xyz", structures, write_results=False)

# Structural optimization and again calculating energy
for s in structures:
    s.calc = calculator
    ucf = UnitCellFilter(s)
    opt = BFGS(ucf)
    opt.run(fmax=0.001)
    s.info["MACE_energy"] = s.get_potential_energy()
    print(s.info["MACE_energy"])

write("vasp/optimized.xyz", structures, write_results=False)
