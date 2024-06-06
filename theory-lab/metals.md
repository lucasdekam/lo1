---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Metals

Most electrodes are made of metal. The electrodes in the fuel cell of a hydrogen car are made of platinum, for example. Here we calculate the binding energy of oxygen to metal electrodes. To do this, we need to calculate the ground state energy of a system.

The metal electrodes in experiments are several millimeters or centimeters big. An atom is less than a nanometer wide. Calculating the density of electrons for the millions of atoms in an electrode obviously takes too long on a computer. Luckily, metal is organized in a *lattice*, and looks similar everywhere.

For this reason, one uses 'periodic boundary conditions': you model only a few atoms, called a 'supercell', and repeat the cell in each direction. This means that electrons on the right of the cell interact with electrons on the left of the next cell. In the picture below, you see such a supercell. 

```{code-cell}
:tags: [hide-input]

from ase.build import bulk
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

metal = bulk('Pt', crystalstructure='fcc', a=4.02)

fig = plt.figure(figsize=(3, 3))
ax = fig.add_subplot()
plot_atoms(metal.repeat((5, 5, 3)), ax)
plt.axis('off')
plt.show()
```


## Density functional theory for periodic systems

Our goal is now to solve the Kohn-Sham equations for the orbitals $\psi_i$ and their energies $\varepsilon_i$ to obtain the ground state energy of a metal. It turns out that it is very efficient to make use of the periodicity of the simulation cell. 

Due to the periodicity, electrons can be expressed as waves that travel through the material. The waves have a frequency in space that is expressed by the three-dimensional *wave vector* $\mathbf{k}=(k_x, k_y, k_z)$. Each frequency is associated with an energy. A plot of the energy for each $\mathbf{k}$ is called the *band structure*. As metals have many electrons, there are also many bands.

```{figure} ../images/bandstructure-bare.png
---
height: 800
name: bandstructure
---
Bandstructure of a platinum fcc(111) surface
```

Expressing electrons in a periodic system as waves is closely related to the techniques of the [Fourier transform](https://youtu.be/spUNpyF58BY?si=Z2UN96n6tIFOa1WV) and [Fourier series](https://youtu.be/r6sGWTCMz2k?si=hFfTSdh6r_OQolRd), which allows any function to be expressed as a sum of waves with various frequencies.

Near a nucleus, the electron wave function changes very fast. Describing such fast changes requires a high spatial frequency $|\mathbf{k}|$. Simulating more $\mathbf{k}$ increases computation time. To speed up calculations, people have defined *pseudopotentials* for many elements. Pseudopotentials include nucleus as well as the core electrons, which shield the charge of the nuclei. This smoothes the wave function, so that we do not need to calculate as many $\mathbf{k}$ points.

```{figure} ../images/pp.png
---
height: 300
name: pseudopotentials
---
Pseudopotential and the corresponding pseudo-wavefunction (both red) comapred to the real potential and real wavefunction. 
```


## Task 1: bulk metal

Let's get started with simulations. We'll use the software [VASP](https://vasp.at/) (Vienna Ab-initio Simulation Package), which runs on the computation cluster of the university. We will write several files that contain the simulation settings (*input files*). You can then upload them to the following website:

[lab.tc.lic.leidenuniv.nl/LGI](https://lab.tc.lic.leidenuniv.nl/LGI/)

VASP will then do the calculation for you.

The first file we need is the `INCAR` file. Make an empty file and name it `INCAR` (not `INCAR.txt`, remove the extension). Enter the following contents:

```
The first line can be used for a comment. For example: "bulk fcc metal"
ENCUT = 300
ISMEAR = 0
SIGMA = 0.05
EDIFF = 1e-6
GGA = RP
```

These are different settings for the algorithm that solves the Kohn-Sham equations. `ENCUT` specifies the maximum frequency of the electron waves used, in terms of an energy through the relation $E=\hbar |\mathbf{k}|$. `EDIFF` is the *convergence criterion* for the self-consistent solving cycles. If the energy of subsequent solutions changes by less than `EDIFF`, the algorithm stops and returns the final result. `ISMEAR` and `SIGMA` spread out the occupation of bands by electrons; this helps convergence. `GGA` sets the exchange-correlation functional to the RPBE functional. 

Next, we need to specify the simulation cell. This is done by making a file `POSCAR` with the following contents:

```
Au
 4.20
     0.5 0.5 0.0
     0.0 0.5 0.5
     0.5 0.0 0.5
 Au
   1
Cartesian
  0 0 0
```

The first line is the chemical symbol of the material. The next line is a scaling factor for the whole cell. The three lines below are normalized vectors of the cell. Here, I chose a cell that contains only one atom, the 'primitive unit cell'. It looks like this:

```{code-cell}
:tags: [hide-input]

from ase.build import bulk
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

metal = bulk('Au', crystalstructure='fcc', a=4.20)

fig = plt.figure(figsize=(3, 3))
ax = fig.add_subplot()
plot_atoms(metal, ax)
plt.axis('off')
plt.show()
```

Next, the `POSCAR` describes how many atoms there are of each type (here: Au, 1). The last lines of the `POSCAR` file set a shift of the atom in the cell. We don't need to shift it, so we set it to zero.

Next, we need a `KPOINTS` file. We cannot calculate all k-points in the bandstructure, so we approximate it by a grid with a number of points. Here's a possible k-points file:

```
K-Points
 0
Gamma
 11 11 11
 0  0  0
```

You can use the documentation on the [VASP wiki](https://www.vasp.at/wiki/index.php/KPOINTS#Regular_k-point_mesh) to understand the different lines.  

Lastly, we need a `POTCAR` file, that contains the pseudopotential that VASP will use. These you will get from the supervisor.

## Task 2: oxygen molecule

A bulk metal is periodic, but an oxygen molecule is not. However, we can still simulate an oxygen molecule by making a large enough simulation cell, so that the molecule in one cell does not interact with the molecule in the next cell. Let's do such a simulation.

INCAR:
```
Oxygen molecule
ENCUT = 400
ISMEAR = 0
SIGMA = 0.05
EDIFF = 1e-06
GGA = RP

EDIFFG = -1e-2
IBRION = 2
ISIF = 2
NSW = 100

ISPIN = 2
```

POSCAR:
```
O
 1.0000000000000000
    10.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   10.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   10.0000000000000000
 O
   2
Selective dynamics
Cartesian
  0  0  0   
  0  0  1   
```

KPOINTS:

The first group of tags is the same as before. We just need a larger `ENCUT` to describe oxygen. The second group is related to geometry relaxation: `IBRION = 2` specifies that we optimize the atoms until the forces on them are below `EDIFFG = -1e-2`, in eV/Å (1 Å, Ångstrom, is 0.1 nm). The last line indicates that we take into account spin. This is because the ground state of the oxygen molecule has two unpaired electrons.



Simulating bulk FCC metal.

Simulating an oxygen molecule. Spin.

Simulating an empty slab. ASE.  (write down lattice constants)

Simulating a slab with adsorbate. (write down lattice constants)

Gibbs free energy of adsorption. Entropy and zero-point energy (a consequence of the Heisenberg uncertainty principle).

Varying simulation parameters. Various surfaces.