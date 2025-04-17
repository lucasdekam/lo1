"""
This script generates VASP input files for a set of configurations.
"""

import os
import argparse
from os import environ
from glob import glob
import numpy as np
from ase.calculators.vasp import Vasp
from ase.io import read
import yaml

KSPACING = 0.2  # per Angstrom


def main(index, output_dir):
    # Set pseudopotential path for ase
    environ["VASP_PP_PATH"] = environ["HOME"] + "/vasp/pps"

    # Read configs
    config_names = glob("ions/frame*.xyz")
    config_names.sort()

    config_names = config_names[:10]  # Limit to first 100 configurations

    # Load DFT parameters
    with open("params.yaml", "r", encoding="utf-8") as f:
        params = yaml.safe_load(f)

    atoms = read("optimized.xyz", index=index)
    directory_name = os.path.join(output_dir)
    os.makedirs(directory_name, exist_ok=True)

    kpts = [1, 1, 1]
    gamma = True
    ispin = 1
    if atoms.info["type"] == "material":
        kpts = [
            max(1, round(2 * np.pi / a / KSPACING))
            for a in atoms.get_cell_lengths_and_angles()[:3]
        ]
    if "spin_polarized" in atoms.info:
        if atoms.info["spin_polarized"]:
            ispin = 2

    calc = Vasp(
        directory=directory_name,
        kpts=kpts,
        gamma=gamma,
        ispin=ispin,
        ncore=4,
        **params[atoms.info["type"]],
    )
    print(f"Writing to directory {directory_name}...", flush=True)
    calc.write_input(atoms)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate VASP input files for a specific configuration."
    )
    parser.add_argument(
        "--index",
        type=int,
        required=True,
        help="Index of the configuration to process.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory for VASP input files.",
    )
    args = parser.parse_args()

    main(args.index, args.output_dir)
