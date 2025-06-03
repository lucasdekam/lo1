# Extra: Quantum chemistry

The following text is inspired by the lecture notes of Christopher Stein (TU Munich) and Katharina Doblhoff-Dier (Leiden Univ.), and the Wikipedia page on [density functional theory](https://en.wikipedia.org/wiki/Density_functional_theory). The introductory article of {cite:t}`BW13` might also be interesting to the reader.

In this practical course we consider electrocatalysis at the scale of atoms. In particular, we are interested in the chemical bond between a reaction intermediate and a catalyst. A chemical bond forms due to the interaction between the electrons of the reaction intermediate and the catalyst. So: to understand chemical bonds, we need to understand electrons.

## The Schrödinger equation

Electrons are very small and very light, much lighter even than atomic nuclei. Electrons do not obey Newton's laws of classical mechanics. Instead, we need to consider quantum mechanics. The central quantity in quantum mechanics is the wave function. If we know the wave function of something, we can use it to calculate any observable quantity.

The wave function of all electrons in a certain system (for example, a molecule) can be found by solving the Schrödinger equation,

$$ \hat{H}(\mathbf{R}) \Psi(\mathbf{R}) = E \Psi(\mathbf{R}) $$

where $\mathbf{R}$ denotes all coordinates of all electrons in the system, $\hat{H}$ is the Hamiltonian of the system, $E$ is the energy and $\Psi$ is the wavefunction. Mathematically, $\Psi(\mathbf{R})$ is an *eigenfunction* of $\hat{H}$, and $E$ the corresponding *eigenvalue*.

In general, we want to calculate the lowest possible energy, the *ground state*. Although it is in principle possible to calculate higher, excited states, these usually do not occur at room temperature.

```{note}
:class: dropdown
This is the time-independent Schrödinger equation. We don't need to consider dynamics of electrons. We also assume that all nuclei are fixed in space. We can make this assumption because the nuclei are much heavier than the electrons. This approximation is also called the Born-Oppenheimer approximation.
```

The Hamiltonian is an operator. You can think of it as a 'recipe' that tells you how the wave function is related to the energy. For the electrons in a molecule, the Hamiltonian is the sum of several parts:

$$ \hat{H} = \hat{T} + \hat{V}_{ne} + \hat{V}_{ee}.$$

The first term, $\hat{T}$, is the kinetic energy of the electrons; $\hat{V}_{ne}$ the interaction between the electrons and the nuclei; and $\hat{V}_{ee}$ is the interaction between the different electrons. The kinetic energy is no problem; there is a closed-form formula ($\hat{T} = -\hbar^2/2m \times  \nabla^2 \Psi$ for each electron). The interaction between electrons and nuclei is also no problem, because we assumed that the nuclei are fixed in space. However, the electron-electron interaction is very complicated (read about the 'many-body problem' on [Wikipedia](https://en.wikipedia.org/wiki/Many-body_problem)).

To solve the Schrödinger equation for many-electron systems (where many is more than three) on a computer within a reasonable calculation time, people have come up with many different approximations, from [Hartree-Fock theory](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method) to [coupled-cluster theory](https://en.wikipedia.org/wiki/Coupled_cluster). All these methods revolve around calculating the wave function. Some of these methods can be extremely accurate, but it takes a lot of time to do such calculations on a computer. It turns out that there is another approach that achieves decent accuracy in much less computation time: the *density functional theory* (DFT).

## Density functional theory

In density functional theory, the central quantity is the electron density, $\rho(\mathbf{r})$, where $\mathbf{r}$ is the vector of spatial coordinates. Walter Kohn and Pierre Hohenberg showed that $\rho(\mathbf{r})$ also contains all relevant information about the system. Later, Walter Kohn and Lu Jeu Sham showed that the problem of calculating the electronic density can be reduced to solving the Schrödinger equation of noninteracting particles (so no $\hat{V}_{ee}$), which made the method very easy to use, while still being accurate. Walter Kohn was awarded the [Nobel Prize in 1998](https://www.nobelprize.org/prizes/chemistry/1998/summary/).

The method of Kohn and Sham comes down to solving $N$ Schrödinger equations, one for each electron. The electrons move in an effective potential $V_\mathrm{eff}(\mathbf{r})$, which depends on the electron density and the positions of the nuclei. The system of these particular $N$ Schrödinger equations is known as the Kohn-Sham equations; they read

$$ \left( -\frac{\hbar^2}{2m} \nabla^2 + V_\mathrm{eff}(\mathbf{r}) \right) \psi_i(\mathbf{r}) = \varepsilon_i \psi_i(\mathbf{r}). $$

Each electron has its own one-electron wavefunction $\psi_i$, known as an *orbital*. These orbitals have an associated energy $\varepsilon_i$. The density is calculated from the orbitals as

$$ \rho(\mathbf{r}) = \sum_{i=1}^{N} |\psi_i(\mathbf{r})|^2. $$

The effective potential depends on the density as follows:

$$ V_\mathrm{eff}(\mathbf{r}) = V_\mathrm{ext}(\mathbf{r}) + \int \frac{\rho(\mathbf{r}')}{4\pi\varepsilon_0 |\mathbf{r} - \mathbf{r}'|} \mathrm d \mathbf{r}' + \frac{\delta E_\mathrm{xc}[\rho(\mathbf{r})]}{\delta \rho(\mathbf{r})}.$$

For this course, it is not important to understand what all these terms mean in detail. In simple words:

* $V_\mathrm{ext}$ (ext for external) describes the electrostatic interaction of the electron with the various nuclei in the system.

* The integral of $\rho$ describes the electrostatic interaction of the electron with the negatively charged electron cloud, consisting of all electrons in the system.

* The *exchange-correlation functional* $E_\mathrm{xc}$ contains corrections to the approximations made by Kohn and Sham. Finding better exchange-correlation functionals is an ongoing topic in research. Different systems might be described better by different exchange-correlation functionals. PBE and r2SCAN are examples of exchange-correlation functionals.

```{note}
:class: dropdown
In case you are wondering what the symbols in the last term mean, here's an explanation. A *function* $f(x)$ takes a number as input, and outputs a number. For example: $f(x)=x^2$ turns $x=2$ into $f(2)=4$. A *functional* $F[f(x)]$ takes an entire function as input, and outputs a number. Functionals also have their own kind of derivative, denoted with a $\delta$.
```

It follows that the Kohn-Sham equations depend on $\rho$ through $V_\mathrm{eff}$. To write down and solve the Kohn-Sham equations, we therefore need to know $\rho$, but to know $\rho$ we need to know the $\psi_i$, which we get by solving the Kohn-Sham equations. How do we solve this circular problem? Luckily, it turns out that when we solve the Kohn-Sham equations for a certain $\rho$, we get out $\psi_i$ that give a 'better' $\rho$. 'Better' means: lower energy, because we are looking for the lowest energy state (ground state). By starting out with a certain guess for $\rho$, we can keep solving the Kohn-Sham equations with an improved $\rho$ every time. When our energy does not change anymore, we know that we are at a minimum: we are *converged*. Such a 'circular' procedure is called a *self-consistent procedure*.

<!-- Now the question remains: how do we solve the Kohn-Sham equations numerically, and what are the $\psi_i$? Depending on whether you have to simulate metal electrodes or porphyrin molecules, you will use different representations of $\psi_i$. These representations will be discussed in the respective pages, and you will use software to use the Kohn-Sham equations self-consistently yourself. -->

<!-- ## Energies

From a DFT calculation we get the ground state energy $E$ of all electrons, which is the sum of energies of all occupied orbitals. You can think of it as the 'bonding energy'. However, the movement of the nuclei also adds energy to the system.

In electrochemical experiments, we usually have a constant temperature $T$ and pressure $P$. The reaction energies are therefore described by Gibbs free energies:

$$G = U - TS + PV$$

with $U$ the internal energy, $S$ the entropy and $V$ the system volume. The internal energy and the entropy are affected by the motion of the nuclei: rotational, translational, and vibrational. The $PV$ term is usually only relevant for gases.

In this practical, we will mostly worry about the electronic energies. The thermodynamic corrections to find the Gibbs energy will be provided.

```{note}
:class: dropdown
Although they are much heavier than electrons, nuclei also experience quantum mechanical effects. By the Heisenberg uncertainty principle, $\Delta x \Delta p \geq \hbar/2$, we cannot precisely determine the position or momentum ('speed') of a nucleus. So even if we say the nucleus is precisely at the lowest energy position, it could have some momentum, and therefore kinetic energy. We could also say it is not moving, but then it might be in a higher energy position.

This 'uncertainty energy' is approximated by the description of a quantum harmonic oscillator, which has energies

$$E_n = (n + \frac12) \hbar \omega.$$

At higher temperatures, higher modes can be occupied.

From a DFT calculation we can calculate the forces on each atom and thereby the vibrational frequencies of the systems. Each vibrational mode has its own $\omega_i$.

The lowest possible energy the harmonic oscillator can have is $\frac12 \hbar \omega$, due to Heisenberg's uncertainty principle. Because this lowest energy is there even at zero Kelvin, it is called the zero-point energy (ZPE).
``` -->
