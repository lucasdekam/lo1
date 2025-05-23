"""
Collect VASP results and write to a single XYZ file.
"""

from os import environ
from ase.calculators.vasp import Vasp
from ase.io import read, write

# Set pseudopotential path for ase
# environ["VASP_PP_PATH"] = environ["HOME"] + "/vasp/pps"
# environ["VASP_PP_PATH"] = "C:\Users\kamlbtde\LocalFiles"

# Read configs
configs = read("optimized.xyz", format="extxyz", index=":")
configs_towrite = []

# Read DFT calculation outputs
for i, atoms in enumerate(configs):
    directory_name = f"configs/{i:05.0f}"
    try:
        calc = Vasp(
            directory=directory_name,
            restart=True,
        )
        atoms.info["DFT_energy"] = calc.get_potential_energy()
        atoms.positions = calc.get_atoms().positions
        atoms.cell = calc.get_atoms().cell
        configs_towrite.append(atoms)
        print(f"Read config {i:05.0f}, " f"energy {atoms.info['DFT_energy']:.3f} eV")
    except RuntimeError:
        print(
            f"Did not read config {i:05.0f}, "
            f"calculation might not have finished. Check manually."
        )

write(
    "dft.xyz",
    configs_towrite,
    format="extxyz",
    columns=[
        "symbols",
        "positions",
    ],
    write_results=False,
)
