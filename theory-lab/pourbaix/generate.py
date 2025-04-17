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
names = [
    "Li",
    "Na",
    "Li2O",
    "Na2O",
    "Co_heterogenite_2h",
    "Co_heterogenite_3r",
    "Co_hydroxide",
    "CoO",
    "CoO2",
    "Co3O4",
]
cod_ids = [
    9008473,
    9008545,
    1010064,
    1010876,
    9009449,
    9009884,
    1010267,
    1528838,
    3000496,
    1538531,
]
materials = [read(name + ".cif", format="cif") for name in names]

for i, m in enumerate(materials):
    m.info["type"] = "material"
    m.info["name"] = names[i]
    m.info["cod_id"] = cod_ids[i]
    if "occupancy" in m.info:
        del m.info["occupancy"]
    if "spacegroup_kinds" in m.arrays:
        del m.arrays["spacegroup_kinds"]

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
    opt.run(fmax=0.02)
    s.info["MACE_energy"] = s.get_potential_energy()
    print(s.info["MACE_energy"])

write("vasp/optimized.xyz", structures, write_results=False)
