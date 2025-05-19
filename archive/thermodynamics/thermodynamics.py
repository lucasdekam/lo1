# https://github.com/patonlab/GoodVibes/blob/master/goodvibes/thermo.py
# https://gaussian.com/wp-content/uploads/dl/thermo.pdf
import numpy as np
import constants as C


def diatomic_gibbs_correction(
    atomic_mass,
    bond_length,
    omega,
    temperature,
    pressure,
):
    """
    Calculate the Gibbs free energy (eV) per particle of a diatomic molecule in the
    ideal gas approximation

    G = E - TS + PV = E - TS + NkT
    """
    internal_energy = (
        translational_energy(temperature)
        + linear_rotational_energy(temperature)
        + vibrational_energy(omega, temperature)
    )
    entropy = (
        translational_entropy(2 * atomic_mass, temperature, pressure)
        + diatomic_rotational_entropy(atomic_mass, bond_length, temperature)
        + vibrational_entropy(omega, temperature)
    )
    gibbs = (
        internal_energy
        - temperature * entropy
        + C.BOLTZMANN_CONSTANT * temperature / C.ELEMENTARY_CHARGE
    )
    return gibbs


def adsorbate_gibbs_correction(
    omega,
    temperature,
):
    """
    Calculate the Gibbs free energy per particle of an adsorbate

    G = E - TS (+ PV)
    PV contribution for the slab and adsorbate is neglected (why?)
    """
    internal_energy = vibrational_energy(omega, temperature)
    entropy = vibrational_entropy(omega, temperature)
    gibbs = internal_energy - temperature * entropy
    return gibbs


def translational_energy(temperature):
    """
    Calculates the translational energy (eV) of an ideal gas molecule,
    i.e., for non-interacting molecules.
    Etrans = 3/2 kT

    Parameter:
    temperature (float): in K.
    """
    energy = 1.5 * C.BOLTZMANN_CONSTANT * temperature
    return energy / C.ELEMENTARY_CHARGE


def linear_rotational_energy(temperature):
    """
    Calculates the rotational energy of a linear molecule (eV)
    Etrans = kT (linear)

    Parameters:
    temperature (float) in K
    """
    energy = C.BOLTZMANN_CONSTANT * temperature
    return energy / C.ELEMENTARY_CHARGE


def vibrational_energy(omega, temperature):
    """
    Vibrational energy evaluation.

    Calculates the vibrational energy contribution (eV).
    Includes ZPE (0K) and thermal contributions.
    Evib = k * Sum(0.5 hw/k + (hw/k)/(e^(hw/kT)-1))

    Parameters:
    omega (float | list): list of angular vibration frequencies
    temperature (float): temperature for calculations to be performed at.

    Returns:
    float: vibrational energy summed up over all modes.
    """
    omega = np.atleast_1d(omega)
    vibrational_temperatures = [
        (C.RED_PLANCK_CONSTANT * w / C.BOLTZMANN_CONSTANT) for w in omega
    ]

    # Error occurs if T is too low when performing np.exp
    for entry in vibrational_temperatures:
        if entry / temperature > np.log(np.finfo(np.float32).max):
            energy = [
                C.BOLTZMANN_CONSTANT * entry * 0.5 for entry in vibrational_temperatures
            ]

    energy = [
        C.BOLTZMANN_CONSTANT * entry * (0.5 + (1 / (np.exp(entry / temperature) - 1)))
        for entry in vibrational_temperatures
    ]

    return np.sum(energy) / C.ELEMENTARY_CHARGE


def zeropoint_energy(omega):
    """
    Vibrational Zero point energy evaluation.

    Calculates the vibrational ZPE (eV)
    EZPE = Sum(0.5 hw/k)

    Parameters:
    omega (list): list of angular vibration frequencies
    temperature (float): temperature for calculations to be performed at.

    Returns:
    float: zero-point vibrational energy summed up over all modes.
    """
    omega = np.atleast_1d(omega)
    energy = [0.5 * C.RED_PLANCK_CONSTANT * w for w in omega]
    return np.sum(energy) / C.ELEMENTARY_CHARGE


def translational_entropy(molecular_mass, temperature, pressure):
    """
    Translational entropy evaluation.

    Calculates the translational entropic contribution (eV/K) of an ideal gas.
    Needs the molecular mass. Convert mass in amu to kg
    Strans = k(ln(2pimkT/h^2)^3/2(kT/P)) + 1 + 3/2)

    Parameters:
    molecular_mass (float): total molecular mass (u) of chemical system.
    conc (float): concentration to perform calculations at.
    temperature (float): temperature for calculations to be performed at.
    pressure (float): pressure of the ideal gas.

    Returns:
    float: translational entropy of chemical system.
    """
    partition = (
        (
            2
            * np.pi
            * molecular_mass
            * C.AMU_TO_KG
            * C.BOLTZMANN_CONSTANT
            * temperature
            / C.PLANCK_CONSTANT**2
        )
        ** (3 / 2)
        * C.BOLTZMANN_CONSTANT
        * temperature
        / pressure
    )
    entropy = C.BOLTZMANN_CONSTANT * (5 / 2 + np.log(partition))
    return entropy / C.ELEMENTARY_CHARGE


def electronic_entropy(multiplicity):
    """
    Calculates the electronic entropic contribution (eV/K) of the molecule
    Selec = k (ln(multiplicity)

    multiplicity (int): spin multiplicity of the molecule
    """
    entropy = C.BOLTZMANN_CONSTANT * np.log(multiplicity)
    return entropy / C.ELEMENTARY_CHARGE


def diatomic_rotational_entropy(atomic_mass, bond_length, temperature):
    """
    Calculates the rotational entropy (eV/K)
    Strans = 0 (atomic) ; R(Ln(q)+1) (linear); R(Ln(q)+3/2) (non-linear)

    Parameters:
    atomic_mass (float): mass (u) of one of the diatomic nuclei (assuming they are the same)
    bond_length (float): bond length (Angstrom) between the two nuclei
    temperature (float): temperature for calculations to be performed at.

    Returns:
    float: rotational entropy of chemical system.
    """
    moment_of_inertia = (
        1 / 2 * atomic_mass * C.AMU_TO_KG * (bond_length * C.ANGSTROM_TO_M) ** 2
    )
    rotational_temp = C.PLANCK_CONSTANT**2 / (
        8 * np.pi**2 * moment_of_inertia * C.BOLTZMANN_CONSTANT
    )
    partition = 1 / 2 * (temperature / rotational_temp)
    entropy = C.BOLTZMANN_CONSTANT * (np.log(partition) + 1)
    return entropy / C.ELEMENTARY_CHARGE


def vibrational_entropy(omega, temperature):
    """
    Entropic contributions (J/K) according to a rigid-rotor
    harmonic-oscillator description for a list of vibrational modes
    Sv = k Sum(hw/(kT(e^(hw/kT)-1) - ln(1-e^(-hw/kT)))

    Parameters:
    omega (list): list of angular frequencies parsed from file.
    temperature (float): temperature for calculations to be performed at.
    freq_scale_factor (float): frequency scaling factor based on level of theory and basis set used.
    fract_modelsys (list): MM frequency scale factors obtained from ONIOM calculations.
    """
    omega = np.atleast_1d(omega)

    exponents = [
        (C.RED_PLANCK_CONSTANT * w) / (C.BOLTZMANN_CONSTANT * temperature)
        for w in omega
    ]
    entropy = [
        C.BOLTZMANN_CONSTANT * entry / (np.exp(entry) - 1)
        - C.BOLTZMANN_CONSTANT * np.log(1 - np.exp(-entry))
        for entry in exponents
    ]
    return np.sum(entropy) / C.ELEMENTARY_CHARGE
